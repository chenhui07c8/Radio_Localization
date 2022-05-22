
% SPP
% get the derivative of the equivalent array response
% steering vector angles are chosen as: PhiBM, PhiRM
function c = get_D_sv_Aeq_SPP(c)
    beamsplit_coe = c.beamsplit_coe;

    % AeqBM
    D_phiBM_AeqBM = zeros(c.NB, c.NM, c.K);
    D_thetaBM_AeqBM = zeros(c.NB, c.NM, c.K);
    D_tauBM_AeqBM = zeros(c.NB, c.NM, c.K); % zeros in SPP
    D_alpha_AeqBM = zeros(c.NB, c.NM, c.K); % zeros in SPP
    D_beta_AeqBM = zeros(c.NB, c.NM, c.K);  % zeros in SPP
    D_gamma_AeqBM = zeros(c.NB, c.NM, c.K); % zeros in SPP
    
    % AeqMB = zeros(c.NM, c.K);
    D_phiBM_AeqMB = zeros(c.NM, c.NB, c.K);
    D_thetaBM_AeqMB = zeros(c.NM, c.NB, c.K);
    D_tauBM_AeqMB = zeros(c.NM, c.NB, c.K); % zeros in SPP
    D_alpha_AeqMB = zeros(c.NM, c.NB, c.K); % 
    D_beta_AeqMB = zeros(c.NM, c.NB, c.K);  % 
    D_gamma_AeqMB = zeros(c.NM, c.NB, c.K); % 
    
    % BM channels
    [D_phiBM_tBM, D_thetaBM_tBM] = get_D_Phi_t(c.phiBM, c.thetaBM);

    for nb = 1:c.NB     % number of antennas
        % BF direction vector
        tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
        
        % BM
        OmegaPsi = c.Bae0'*c.RotB'*c.tBM - (c.Bae0'*tBFtemp)*beamsplit_coe;
        AeqMat = 1/sqrt(c.DsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
%         AeqBM(nb, :) = abs(sum(AeqMat,1));  % imaginary part is 0 after summation
        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_phiBM_tBM));
        D_phiBM_AeqBM(nb, :, :) = repmat(sum(D_AeqMat,1), c.NM, 1);
        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_thetaBM_tBM));
        D_thetaBM_AeqBM(nb, :, :) = repmat(sum(D_AeqMat,1), c.NM, 1);
    end
    % local DOD/DOA: global = Rotm*local; local = Rotm^-1*global
    
    % MB channels
    OM = c.OM;
    [D_alpha_Rot, D_beta_Rot, D_gamma_Rot] = get_D_Alpha_Rot(OM);
    [D_phiBM_tBM, D_thetaBM_tBM] = get_D_Phi_t(c.phiBM, c.thetaBM);
    D_phiBM_tMB = -D_phiBM_tBM;
    D_thetaBM_tMB = -D_thetaBM_tBM;
    
    for nm = 1:c.NM     % number of antennas
        % get beamforming gain
        tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.

        % MB
        OmegaPsi = c.Mae0'*c.RotM'*c.tMB - (c.Mae0'*tBFtemp)*beamsplit_coe;
        AeqMat = 1/sqrt(c.DsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
%         AeqMB(nm, :) = abs(sum(AeqMat,1));

        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_alpha_Rot'*c.tMB));
        D_alpha_AeqMB(nm, :, :) = repmat(sum(D_AeqMat,1), c.NB, 1);   % sum over Nsa antenna elements.

        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_beta_Rot'*c.tMB));
        D_beta_AeqMB(nm, :, :) = repmat(sum(D_AeqMat,1), c.NB, 1);

        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_gamma_Rot'*c.tMB));
        D_gamma_AeqMB(nm, :, :) = repmat(sum(D_AeqMat,1), c.NB, 1);

        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*c.RotM'*D_phiBM_tMB));
        D_phiBM_AeqMB(nm, :, :) = repmat(sum(D_AeqMat,1), c.NB, 1);

        D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*c.RotM'*D_thetaBM_tMB));
        D_thetaBM_AeqMB(nm, :, :) = repmat(sum(D_AeqMat,1), c.NB, 1);
    end
    
%     c.AeqBM = AeqBM;
    c.D_phiBM_AeqBM = D_phiBM_AeqBM;
    c.D_thetaBM_AeqBM = D_thetaBM_AeqBM;
    c.D_tauBM_AeqBM = D_tauBM_AeqBM;
    c.D_alpha_AeqBM = D_alpha_AeqBM;
    c.D_beta_AeqBM = D_beta_AeqBM;
    c.D_gamma_AeqBM = D_gamma_AeqBM;

%     c.AeqBR = AeqBR;
    
%     c.AeqMB = AeqMB;
    c.D_phiBM_AeqMB = D_phiBM_AeqMB;
    c.D_thetaBM_AeqMB = D_thetaBM_AeqMB;
    c.D_tauBM_AeqMB = D_tauBM_AeqMB;
    c.D_alpha_AeqMB = D_alpha_AeqMB;
    c.D_beta_AeqMB = D_beta_AeqMB;
    c.D_gamma_AeqMB = D_gamma_AeqMB;

    % AeqMR = zeros(c.NM, c.K);
    D_phiRM_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
    D_thetaRM_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
    D_tauRM_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K); % zeros in SPP
    D_alpha_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K); % zeros in SPP
    D_beta_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);  % zeros in SPP
    D_gamma_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K); % zeros in SPP

    D_phiRM_AeqMR = zeros(c.NM, c.NR, c.K);
    D_thetaRM_AeqMR = zeros(c.NM, c.NR, c.K);
    D_tauRM_AeqMR = zeros(c.NM, c.NR, c.K); % zeros in SPP
    D_alpha_AeqMR = zeros(c.NM, c.NR, c.K);
    D_beta_AeqMR = zeros(c.NM, c.NR, c.K);
    D_gamma_AeqMR = zeros(c.NM, c.NR, c.K);
    if (c.LR > 0 && c.NsaR > 1)
        % AeqMR = zeros(c.NM, c.K);
        
        [D_phiRM_tRM, D_thetaRM_tRM] = get_D_Phi_t(c.phiRM, c.thetaRM);
        D_phiRM_tMR = -D_phiRM_tRM;
        D_thetaRM_tMR = -D_thetaRM_tRM;
        
        % RM channel (AeqBRM)
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
                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*c.RotM'*D_phiRM_tRM));
                    D_phiRM_AeqBRM(nb, nr, nm, :) = sum(D_AeqMat,1);

                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*c.RotM'*D_thetaRM_tRM));
                    D_thetaRM_AeqBRM(nb, nr, nm, :) = sum(D_AeqMat,1);
                end
            end
        end

        
        % MR channel (AeqMR)
        for nm = 1:c.NM     % number of antennas
            % get beamforming gain
            tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.

            % MR
            OmegaPsi = c.Mae0'*c.RotM'*c.tMR - (c.Mae0'*tBFtemp)*beamsplit_coe;
            
            AeqMat = 1/sqrt(c.DsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_alpha_Rot'*c.tMR));
            D_alpha_AeqMR(nm, :, :) = repmat(sum(D_AeqMat,1), c.NR, 1);   % sum over Nsa antenna elements.

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_beta_Rot'*c.tMR));
            D_beta_AeqMR(nm, :, :) = repmat(sum(D_AeqMat,1), c.NR, 1);

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_gamma_Rot'*c.tMR));
            D_gamma_AeqMR(nm, :, :) = repmat(sum(D_AeqMat,1), c.NR, 1);

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*c.RotM'*D_phiRM_tMR));
            D_phiRM_AeqMR(nm, :, :) = repmat(sum(D_AeqMat,1), c.NR, 1);

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*c.RotM'*D_thetaRM_tMR));
            D_thetaRM_AeqMR(nm, :, :) = repmat(sum(D_AeqMat,1), c.NR, 1);
        end
    end
            
    c.D_phiRM_AeqBRM = D_phiRM_AeqBRM;
    c.D_thetaRM_AeqBRM = D_thetaRM_AeqBRM;
    c.D_tauRM_AeqBRM = D_tauRM_AeqBRM;
    c.D_alpha_AeqBRM = D_alpha_AeqBRM;
    c.D_beta_AeqBRM = D_beta_AeqBRM;
    c.D_gamma_AeqBRM = D_gamma_AeqBRM;

%     c.AeqMR = AeqMR;
    c.D_phiRM_AeqMR = D_phiRM_AeqMR;
    c.D_thetaRM_AeqMR = D_thetaRM_AeqMR;
    c.D_tauRM_AeqMR = D_tauRM_AeqMR;
    c.D_alpha_AeqMR = D_alpha_AeqMR;
    c.D_beta_AeqMR = D_beta_AeqMR;
    c.D_gamma_AeqMR = D_gamma_AeqMR;
        
    
    % NLOS has only SPP model (each reflector is an independent object).
    if (c.LC > 0)
        % AeqBC
        D_phiBC_AeqBC = zeros(c.NB, c.K, c.LC);
        D_thetaBC_AeqBC = zeros(c.NB, c.K, c.LC);
        D_tauBC_AeqBC = zeros(c.NB, c.K, c.LC); % zeros in SPP
        D_alpha_AeqBC = zeros(c.NB, c.K, c.LC); % zeros in SPP
        D_beta_AeqBC = zeros(c.NB, c.K, c.LC);  % zeros in SPP
        D_gamma_AeqBC = zeros(c.NB, c.K, c.LC); % zeros in SPP

        % AeqMC = zeros(c.NM, c.K, c.LC);
        D_phiCM_AeqMC = zeros(c.NM, c.K, c.LC);
        D_thetaCM_AeqMC = zeros(c.NM, c.K, c.LC);
        D_tauCM_AeqMC = zeros(c.NM, c.K, c.LC); % zeros in SPP
        D_alpha_AeqMC = zeros(c.NM, c.K, c.LC);
        D_beta_AeqMC = zeros(c.NM, c.K, c.LC);
        D_gamma_AeqMC = zeros(c.NM, c.K, c.LC);
        
        % BC channel
        for l = 1:c.LC
        [D_phiBC_tBC, D_thetaBC_tBC] = get_D_Phi_t(c.phiBC(l), c.thetaBC(l));
            for nb = 1:c.NB     % number of antennas
                % BF direction vector
                tBFtemp = get_dir_from_angle(c.BFmatB(:, nb));  % get the BF direction vector for the b-th SA.
                % BM
                OmegaPsi = c.Bae0'*c.RotB'*c.tBC(:,l) - (c.Bae0'*tBFtemp)*c.beamsplit_coe;
                AeqMat = 1/sqrt(c.DsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_phiBC_tBC));
                D_phiBC_AeqBC(nb, :, l) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_thetaBC_tBC));
                D_thetaBC_AeqBC(nb, :, l) = sum(D_AeqMat,1);
            end
        end
        % MC channel
        OM = c.OM;
        [D_alpha_Rot, D_beta_Rot, D_gamma_Rot] = get_D_Alpha_Rot(OM);
        for l = 1:c.LC
            [D_phiCM_tMC, D_thetaCM_tMC] = get_D_Phi_t(c.phiMC(l), c.thetaMC(l));
            
            for nm = 1:c.NM     % number of antennas
                % get beamforming gain
                tBFtemp = get_dir_from_angle(c.BFmatM(:, nm));  % get the BF direction vector for the b-th SA.

                % MB
                OmegaPsi = c.Mae0'*c.RotM'*c.tMC(:,l) - (c.Mae0'*tBFtemp)*c.beamsplit_coe;
                AeqMat = 1/sqrt(c.DsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
        %         AeqMB(nm, :) = abs(sum(AeqMat,1));

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_alpha_Rot'*c.tMC(:,l)));
                D_alpha_AeqMC(nm, :, l) = sum(D_AeqMat,1);   % sum over Nsa antenna elements.

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_beta_Rot'*c.tMC(:,l)));
                D_beta_AeqMC(nm, :, l) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_gamma_Rot'*c.tMC(:,l)));
                D_gamma_AeqMC(nm, :, l) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*c.RotM'*D_phiCM_tMC));
                D_phiCM_AeqMC(nm, :, l) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*c.RotM'*D_thetaCM_tMC));
                D_thetaCM_AeqMC(nm, :, l) = sum(D_AeqMat,1);
            end
        end

        % BC channel
        c.D_phiBC_AeqBC = D_phiBC_AeqBC;
        c.D_thetaBC_AeqBC = D_thetaBC_AeqBC;
        c.D_tauBC_AeqBC = D_tauBC_AeqBC; % zeros in SPP
        c.D_alpha_AeqBC = D_alpha_AeqBC; % zeros in SPP
        c.D_beta_AeqBC = D_beta_AeqBC;  % zeros in SPP
        c.D_gamma_AeqBC = D_gamma_AeqBC; % zeros in SPP

        % MC channel
        % AeqMC = zeros(c.NM, c.K, c.LC);
        c.D_phiCM_AeqMC = D_phiCM_AeqMC;
        c.D_thetaCM_AeqMC = D_thetaCM_AeqMC;
        c.D_tauCM_AeqMC = D_tauCM_AeqMC; % zeros in SPP
        c.D_alpha_AeqMC = D_alpha_AeqMC;
        c.D_beta_AeqMC = D_beta_AeqMC;
        c.D_gamma_AeqMC = D_gamma_AeqMC;
    end
    

end