% Get the array factor for BS/RIS/MD
function c = get_Aeq_matrix(c)

%************** LOS Channel ***************
    
    % check if BF_angle matches NB, if not, set default as all 0
    if size(c.BFmatB,2) ~= c.NB
        disp('BS beamforming angle dimenion not matching, set all as 0s.')
        c.BFmatB = zeros(2, prod(c.NB_dim))+[0 0]';
    end
    if size(c.BFmatM,2) ~= c.NM
        disp('MD beamforming angle dimenion not matching, set all as 0s.')
        c.BFmatM = zeros(2, prod(c.NM_dim))+[0 0]';
    end
    
    % AOSA structure, SA parameters
    c.Bae0 = c.DaeB*c.dant*get_array_layout(c.NsaB_dim);
    c.Mae0 = c.DaeM*c.dant*get_array_layout(c.NsaM_dim);
    
    AeqBM = zeros(c.NB, c.NM, c.K);
    AeqMB = zeros(c.NM, c.NB ,c.K);
    if c.AoSA_type == "SPP" || c.AoSA_type == "PPP" % S for SA phase, P for SA angle, P for AE beamforming
        % BS beamforming, BM channels
        PsiBM = c.Bae0'*c.RotB'*c.tBM; 
        for nb = 1:c.NB
            % BF direction vector
            tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
            % BM
            OmegaPsi = PsiBM - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
            AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
            AeqBM(nb, :, :) = repmat(real(sum(AeqMat,1)), c.NM, 1);  % imaginary part is 0 after summation
        end
        % local DOD/DOA: global = Rotm*local; local = Rotm^-1*global

        % MD beamforming, MB channels
        PsiMB = c.Mae0'*c.RotM'*c.tMB;
        for nm = 1:c.NM     % number of antennas
            % get beamforming gain
            tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.
            % MB
            OmegaPsi = PsiMB - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
            AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
            AeqMB(nm, :, :) = repmat(real(sum(AeqMat,1)), c.NB, 1);
        end
    elseif c.AoSA_type == "SSP"
        % BS beamforming, BM channels
        for nb = 1:c.NB     % number of antennas
            for nm = 1:c.NM
                tbm = (c.M(:, nm)-c.B(:, nb))/norm(c.M(:, nm)-c.B(:, nb));
                Psibm = c.Bae0'*c.RotB'*tbm;
                % BF direction vector
                tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
                % Psi BM
                OmegaPsi = Psibm - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
                AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                AeqBM(nb, nm, :) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
            end
        end
        % MD beamforming, MB channels
        for nm = 1:c.NM     % number of antennas
            % get beamforming gain
            tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.
            for nb = 1:c.NB
                tmb = (c.B(:, nb)-c.M(:, nm))/norm(c.B(:, nb)-c.M(:, nm));
                Psimb = c.Mae0'*c.RotM'*tmb;
                % MB
                OmegaPsi = Psimb - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
                AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                AeqMB(nm, nb, :) = real(sum(AeqMat,1));
            end
        end
    elseif c.AoSA_type == "SSS"
        AeqBM = ones(size(AeqBM));
        AeqMB = ones(size(AeqMB));
    end
    
    AeqL = zeros(c.NB, c.NM, c.K);
    for k = 1:c.K
        AeqL(:,:,k) = AeqBM(:,:,k).*AeqMB(:,:,k)';
    end
    c.AeqBM = AeqBM;
    c.AeqMB = AeqMB;
    c.AeqL = AeqL;
    
    
%************** RIS channel ***************
    if(c.LR > 0)
        
        % check if BF_angle matches NB, if not, set default as all 0
        if size(c.BFmatRB,2) ~= c.NR
%             disp('RIS beamforming angle dimenion not matching, set all as 0s.')
            c.BFmatRB = zeros(2, prod(c.NR_dim))+[0 0]';
        end
        if size(c.BFmatRM,2) ~= c.NR
%             disp('RIS beamforming angle dimenion not matching, set all as 0s.')
            c.BFmatRM = zeros(2, prod(c.NR_dim))+[0 0]';
        end

        % AOSA BS-RIS-MD RIS channel
        c.Rae0 = c.DsaR*c.dant*get_array_layout(c.NsaR_dim);

        % Effective array response
        AeqBR = zeros(c.NB, c.NR, c.K);
        AeqMR = zeros(c.NM, c.NR, c.K);
        AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
        if c.AoSA_type == "SPP"
            % RIS beamforming, BRM channels
            PsiRB = c.Rae0'*c.RotR'*c.tRB;
            PsiRM = c.Rae0'*c.RotR'*c.tRM;
            for nb = 1:c.NB     % number of antennas
                for nm = 1:c.NM
                    for nr = 1:c.NR       
                        tBFtemp = get_dir_from_angle(c.BFmatRB(:, nr));  % get the BF direction vector for the b-th SA.
                        OmegaPsiRB = PsiRB - (c.Rae0'*tBFtemp)*c.beamsplit_coe;
                        tBFtemp = get_dir_from_angle(c.BFmatRM(:, nr));  % get the BF direction vector for the b-th SA.
                        OmegaPsiRM = PsiRM - (c.Rae0'*tBFtemp)*c.beamsplit_coe;

                        AeqMat = exp(1j*2*pi./c.lambdak.*(OmegaPsiRB + OmegaPsiRM));
                        AeqBRM(nb, nr, nm, :) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
                    end
                end
            end
            
            % BS beamforming, BR channels
            PsiBR = c.Bae0'*c.RotB'*c.tBR;
            for nb = 1:c.NB     % number of antennas
                % BF direction vector
                tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
                % BR
                OmegaPsi = PsiBR - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
                AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                AeqBR(nb, :, :) = repmat(real(sum(AeqMat,1)), c.NR, 1);  % imaginary part is 0 after summation
            end
            % local DOD/DOA: global = Rotm*local; local = Rotm^-1*global

            % MD beamforming, MR channels
            PsiMR = c.Mae0'*c.RotM'*c.tMR; 
            for nm = 1:c.NM     % number of antennas
                % get beamforming gain
                tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.
                % MR
                OmegaPsi = PsiMR - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
                AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                AeqMR(nm, :, :) = repmat(real(sum(AeqMat, 1)), c.NR, 1);
            end
        elseif c.AoSA_type == "SSP"
            
            % BS beamforming, BR channels
            for nb = 1:c.NB     % number of antennas
                for nr = 1:c.NR
                    tbr = (c.R(:, nr)-c.B(:, nb))/norm(c.R(:, nr)-c.B(:, nb));
                    Psibr = c.Bae0'*c.RotB'*tbr;
                    % BF direction vector
                    tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
                    % Psi BR
                    OmegaPsi = Psibr - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
                    AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                    AeqBR(nb, nr, :) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
                end
            end
            
            % MD beamforming, MR channels
            for nm = 1:c.NM     % number of antennas
                % get beamforming gain
                tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.
                for nr = 1:c.NR
                    tmr = (c.R(:, nr)-c.M(:, nm))/norm(c.R(:, nr)-c.M(:, nm));
                    Psimr = c.Mae0'*c.RotM'*tmr;
                    % MB
                    OmegaPsi = Psimr - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
                    AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                    AeqMR(nm, nr, :) = real(sum(AeqMat,1));
                end
            end
            
            % RIS beamforming, BRM channels
            for nb = 1:c.NB     % number of antennas
                for nm = 1:c.NM
                    for nr = 1:c.NR
                        trb = (c.B(:, nb)-c.R(:, nr))/norm(c.B(:, nb)-c.R(:, nr));
                        Psirb = c.Rae0'*c.RotR'*trb;            
                        tBFtemp = get_dir_from_angle(c.BFmatRB(:, nr));  % get the BF direction vector for the b-th SA.
                        OmegaPsirb = Psirb - (c.Rae0'*tBFtemp)*c.beamsplit_coe;

                        trm = (c.M(:, nm)-c.R(:, nr))/norm(c.M(:, nm)-c.R(:, nr));
                        Psirm = c.Rae0'*c.RotR'*trm;            
                        tBFtemp = get_dir_from_angle(c.BFmatRM(:, nr));  % get the BF direction vector for the b-th SA.
                        OmegaPsirm = Psirm - (c.Rae0'*tBFtemp)*c.beamsplit_coe;

                        AeqMat = exp(1j*2*pi./c.lambdak.*(OmegaPsirb + OmegaPsirm));
                        AeqBRM(nb, nr, nm, :) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
                    end
                end
            end
        end
        AeqR = zeros(c.NB, c.NR, c.NM, c.K);
        for nb = 1:c.NB     % number of antennas
            for nm = 1:c.NM
                for nr = 1:c.NR
                    AeqR(nb, nr, nm, :) = reshape(AeqBRM(nb, nr, nm, :), [c.K, 1]).*reshape(AeqBR(nb, nr, :), [c.K 1])...
                        .*reshape(AeqMR(nm, nr, :), [c.K 1]);  % imaginary part is 0 after summation
                end
            end
        end
        c.AeqR = AeqR;
        c.AeqBR = AeqBR;
        c.AeqMR = AeqMR;
        c.AeqBRM = AeqBRM;
    end
    
    
%************** NLOS channel ***************
    if(c.LC > 0)
        AeqBC = zeros(c.NB, c.K, c.LC);
        AeqMC = zeros(c.NM, c.K, c.LC);
        
        if c.AoSA_type == "SPP"
            % BC channels
            PsiBC = c.Bae0'*c.RotB'*c.tBC;
            for nb = 1:c.NB     % number of antennas
                % BF direction vector
                tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.

                for lc = 1:c.LC
                    % BC
                    OmegaPsi = PsiBC(:,lc) - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
                    AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                    AeqBC(nb, :, lc) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
                end
            end

            % MC channels
            PsiMC = c.Mae0'*c.RotM'*c.tMC;
            for nm = 1:c.NM     % number of antennas
                % get beamforming gain
                tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.
                for lc = 1:c.LC
                    % MC
                    OmegaPsi = PsiMC(:,lc) - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
                    AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                    AeqMC(nm, :, lc) = real(sum(AeqMat,1));
                end
            end
        elseif c.AoSA_type == "SSP"
            % BC channels
            for nb = 1:c.NB     % number of antennas
                for lc = 1:c.LC
                    tbc = (c.PC(:, lc)-c.B(:, nb))/norm(c.PC(:, lc)-c.B(:, nb));
                    Psibc = c.Bae0'*c.RotB'*tbc;
                    % BF direction vector
                    tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
                    % Psi BM
                    OmegaPsi = Psibc - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
                    AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                    AeqBC(nb, :, lc) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
                end
            end
            % MC channels
            for nm = 1:c.NM     % number of antennas
                % get beamforming gain
                tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.
                for lc = 1:c.LC
                    tmc = (c.PC(:, lc)-c.M(:, nm))/norm(c.PC(:, lc)-c.M(:, nm));
                    Psimc = c.Mae0'*c.RotM'*tmc;
                    % MC
                    OmegaPsi = Psimc - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
                    AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                    AeqMC(nm, :, lc) = real(sum(AeqMat,1));
                end
            end
        end
        c.AeqBC = AeqBC;
        c.AeqMC = AeqMC;
    end
    
    
    
end