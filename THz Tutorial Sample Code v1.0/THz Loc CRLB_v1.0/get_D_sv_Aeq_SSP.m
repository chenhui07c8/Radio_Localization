
% SSP
% get the derivative of the equivalent array response
% steering vector angles are chosen as: Phibm, Phirm
function c = get_D_sv_Aeq_SSP(c)

    D_RotM_pm = get_rotation_over_position(c);
    B = c.B;
    M = c.M;
    [D_alpha_Rot, D_beta_Rot, D_gamma_Rot] = get_D_Alpha_Rot(c.OM);
    phiBM = c.phiBM;
    thetaBM = c.thetaBM;
    
    
    
    D_phiBM_AeqBM = zeros(c.NB, c.NM, c.K);
    D_thetaBM_AeqBM = zeros(c.NB, c.NM, c.K);
    D_tauBM_AeqBM = zeros(c.NB, c.NM, c.K);
    D_alpha_AeqBM = zeros(c.NB, c.NM, c.K);
    D_beta_AeqBM = zeros(c.NB, c.NM, c.K);
    D_gamma_AeqBM = zeros(c.NB, c.NM, c.K);

    D_phiBM_AeqMB = zeros(c.NM, c.NB, c.K);
    D_thetaBM_AeqMB = zeros(c.NM, c.NB, c.K);
    D_tauBM_AeqMB = zeros(c.NM, c.NB, c.K);
    D_alpha_AeqMB = zeros(c.NM, c.NB, c.K);
    D_beta_AeqMB = zeros(c.NM, c.NB, c.K);
    D_gamma_AeqMB = zeros(c.NM, c.NB, c.K);

    
    beamsplit_coe = c.beamsplit_coe;
    % BM & MB channel: AeqBM, AeqMB
    
    [D_phiBM_tBM, D_thetaBM_tBM] = get_D_Phi_t(c.phiBM, c.thetaBM);
%     D_phiBM_tMB = -D_phiBM_tBM;
%     D_thetaBM_tMB = -D_thetaBM_tBM;

    D_thetaBM_taubm = zeros(c.NM, c.NB);
    D_phiBM_taubm = zeros(c.NM, c.NB);
    D_tauBM_taubm = zeros(c.NM, c.NB);
    for b = 1:c.NB     % number of antennas at BS
        tBFtempB = get_dir_from_angle(c.BFmatB(:, b));  % get the BF direction vector for the b-th SA.
        for m = 1:c.NM
            tBFtempM = get_dir_from_angle(c.BFmatM(:, m));  % get the BF direction vector for the b-th SA.
            % D_svbm_taubm: derivative of sv w.r.t taubm.
            taubm = norm(B(:,b) - M(:,m))/c.lightspeed;
            dbm = norm(B(:,b)-M(:,m));
            % phi
            t_temp = [cosd(thetaBM)*sind(phiBM); -cosd(thetaBM)*cosd(phiBM); 0];
            
            D_phiBM_taubm(m, b) = c.dBM/c.lightspeed/dbm*((M(:,m)-c.PM-B(:,b)+c.PB)'*t_temp);
            % theta
            t_temp = [sind(thetaBM)*cosd(phiBM); sind(thetaBM)*sind(phiBM); -cosd(thetaBM)];
            D_thetaBM_taubm(m, b) = c.dBM/c.lightspeed/dbm*((M(:,m)-c.PM-B(:,b)+c.PB)'*t_temp);
            % tau
            Gbm2 = (M(:,m)-c.PM - B(:,b)+c.PB)'*c.tBM;
            D_tauBM_taubm(m, b) = 1/dbm*(c.dBM - Gbm2);
            D_PhiM_taubm = 1/c.lightspeed/dbm*((M(:,m)-c.PM-B(:,b)+c.PB)'*D_RotM_pm{m} - c.dBM*(c.tBM'*D_RotM_pm{m}));
                
            % ST direction vector
            tbm = (c.M(:, m) - c.B(:, b)) / norm(c.M(:, m) - c.B(:, b));
            
            % BM channel
            OmegaPsi = c.Bae0'*c.RotB'*tbm - (c.Bae0'*tBFtempB)*beamsplit_coe;
            AeqMat = 1/sqrt(c.NsaB)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
    %         AeqBM(nb, :) = abs(sum(AeqMat,1));  % imaginary part is 0 after summation
            D_phiBM_tbm = (D_phiBM_tBM*c.dBM*dbm - (M(:,m)-B(:,b))*c.lightspeed*D_phiBM_taubm(m, b))/dbm^2;
            D_thetaBM_tbm = (D_thetaBM_tBM*c.dBM*dbm - (M(:,m)-B(:,b))*c.lightspeed*D_thetaBM_taubm(m, b))/dbm^2;
            D_tauBM_tbm = (c.tBM*c.lightspeed*dbm - (M(:,m)-B(:,b))*c.lightspeed*D_tauBM_taubm(m, b))/dbm^2;
            D_alpha_tbm = (D_alpha_Rot'*c.M(:,m)*dbm - (M(:,m)-B(:,b))*c.lightspeed*D_PhiM_taubm(1))/dbm^2;
            D_beta_tbm = (D_beta_Rot'*c.M(:,m)*dbm - (M(:,m)-B(:,b))*c.lightspeed*D_PhiM_taubm(2))/dbm^2;
            D_gamma_tbm = (D_gamma_Rot'*c.M(:,m)*dbm - (M(:,m)-B(:,b))*c.lightspeed*D_PhiM_taubm(3))/dbm^2;

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_phiBM_tbm));
            D_phiBM_AeqBM(b, m, :) = sum(D_AeqMat,1);

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_thetaBM_tbm));
            D_thetaBM_AeqBM(b, m, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_tauBM_tbm));
            D_tauBM_AeqBM(b, m, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_alpha_tbm));
            D_alpha_AeqBM(b, m, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_beta_tbm));
            D_beta_AeqBM(b, m, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Bae0'*c.RotB'*D_gamma_tbm));
            D_gamma_AeqBM(b, m, :) = sum(D_AeqMat,1);
            
            % MB channel
            OmegaPsi = c.Mae0'*c.RotM'*c.tMB - (c.Mae0'*tBFtempM)*beamsplit_coe;
            AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
            % Difference: let tmb = -Rm^T*tbm (should be tilde_tmb)
            D_phiBM_tmb = -c.RotM'*D_phiBM_tbm;
            D_thetaBM_tmb = -c.RotM'*D_thetaBM_tbm;
            D_tauBM_tmb = -c.RotM'*D_tauBM_tbm;
            D_alpha_tmb = - D_alpha_Rot'*(tbm) - c.RotM'*D_alpha_tbm;
            D_beta_tmb = - D_beta_Rot'*(tbm) - c.RotM'*D_beta_tbm;
            D_gamma_tmb = - D_gamma_Rot'*(tbm) - c.RotM'*D_gamma_tbm;
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_phiBM_tmb));
            D_phiBM_AeqMB(m, b, :) = sum(D_AeqMat,1);

            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_thetaBM_tmb));
            D_thetaBM_AeqMB(m, b, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_tauBM_tmb));
            D_tauBM_AeqMB(m, b, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_alpha_tmb));
            D_alpha_AeqMB(m, b, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_beta_tmb));
            D_beta_AeqMB(m, b, :) = sum(D_AeqMat,1);
            
            D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_gamma_tmb));
            D_gamma_AeqMB(m, b, :) = sum(D_AeqMat,1);
            
            
        end
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
    
    
% MR channel: AeqMR
    if (c.LR > 0)
        R = c.R;
        phiRM = c.phiRM;
        thetaRM = c.thetaRM;
        D_phiRM_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
        D_thetaRM_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
        D_tauRM_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K); 
        D_alpha_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
        D_beta_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K); 
        D_gamma_AeqBRM = zeros(c.NB, c.NR, c.NM, c.K);
        
        D_phiRM_AeqMR = zeros(c.NM, c.NR, c.K);
        D_thetaRM_AeqMR = zeros(c.NM, c.NR, c.K);
        D_tauRM_AeqMR = zeros(c.NM, c.NR, c.K);
        D_alpha_AeqMR = zeros(c.NM, c.NR, c.K);
        D_beta_AeqMR = zeros(c.NM, c.NR, c.K);
        D_gamma_AeqMR = zeros(c.NM, c.NR, c.K);

        [D_phiRM_tRM, D_thetaRM_tRM] = get_D_Phi_t(c.phiRM, c.thetaRM);
        D_phiRM_tMR = -D_phiRM_tRM;
        D_thetaRM_tMR = -D_thetaRM_tRM;
        D_thetaRM_taurm = zeros(c.NM, c.NR);
        D_phiRM_taurm = zeros(c.NM, c.NR);
        D_tauRM_taurm = zeros(c.NM, c.NR);
        for m = 1:c.NM     % number of antennas at BS
            for r = 1:c.NR
                tBFtempM = get_dir_from_angle(c.BFmatM(:, m));  % get the BF direction vector for the b-th SA.

                % D_svbm_taubm: derivative of sv w.r.t taubm.
                taurm = norm(R(:,r) - M(:,m))/c.lightspeed;
                drm = norm(R(:,r)-M(:,m));
                % theta
                t_temp = [sind(thetaRM)*cosd(phiRM); sind(thetaRM)*sind(phiRM); -cosd(thetaRM)];
                D_thetaRM_taurm(m, r) = c.dRM/c.lightspeed/drm*((M(:,m)-c.PM-R(:,r)+c.PR)'*t_temp);
                % phi
                t_temp = [cosd(thetaRM)*sind(phiRM); -cosd(thetaRM)*cosd(phiRM); 0];
                D_phiRM_taurm(m, r) = c.dRM/c.lightspeed/drm*((M(:,m)-c.PM-R(:,r)+c.PR)'*t_temp);
                % tau
                Gbm2 = (M(:,m)-c.PM - R(:,r)+c.PR)'*c.tRM;
                D_tauRM_taurm(m, r) = 1/drm*(c.dRM - Gbm2);
                D_PhiM_taurm = 1/c.lightspeed/drm*((M(:,m)-c.PM-R(:,r)+c.PR)'*D_RotM_pm{m} - c.dRM*(c.tRM'*D_RotM_pm{m}));

                % ST direction vector
                trm = (c.M(:, m) - c.R(:, r)) / norm(c.M(:, m) - c.R(:, r));
                            
                % MR channel
                D_phiRM_trm = (D_phiRM_tRM*c.dRM*drm - (M(:,m)-R(:,r))*c.lightspeed*D_thetaRM_taurm(m, r))/drm^2;
                D_thetaRM_trm = (D_thetaRM_tRM*c.dRM*drm - (M(:,m)-R(:,r))*c.lightspeed*D_thetaRM_taurm(m, r))/drm^2;
                D_tauRM_trm = (c.tBM*c.lightspeed*drm - (M(:,m)-R(:,r))*c.lightspeed*D_tauRM_taurm(m, r))/drm^2;
                D_alpha_trm = (D_alpha_Rot'*c.M(:,m)*drm - (M(:,m)-R(:,r))*c.lightspeed*D_PhiM_taurm(1))/drm^2;
                D_beta_trm = (D_beta_Rot'*c.M(:,m)*drm - (M(:,m)-R(:,r))*c.lightspeed*D_PhiM_taurm(2))/drm^2;
                D_gamma_trm = (D_gamma_Rot'*c.M(:,m)*drm - (M(:,m)-R(:,r))*c.lightspeed*D_PhiM_taurm(3))/drm^2;

                OmegaPsi = c.Mae0'*c.RotM'*c.tMR - (c.Mae0'*tBFtempM)*beamsplit_coe;
                AeqMat = 1/sqrt(c.NsaM)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
                % Difference: let tmr = -Rm^T*trm (should be tilde_tmr)
                D_phiRM_tmr = -c.RotM'*D_phiRM_trm;
                D_thetaRM_tmr = -c.RotM'*D_thetaRM_trm;
                D_tauRM_tmr = -c.RotM'*D_tauRM_trm;
                D_alpha_tmr = - D_alpha_Rot'*(trm) - c.RotM'*D_alpha_trm;
                D_beta_tmr = - D_beta_Rot'*(trm) - c.RotM'*D_beta_trm;
                D_gamma_tmr = - D_gamma_Rot'*(trm) - c.RotM'*D_gamma_trm;

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_phiRM_tmr));
                D_phiRM_AeqMR(m, r, :) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_thetaRM_tmr));
                D_thetaRM_AeqMR(m, r, :) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_tauRM_tmr));
                D_tauRM_AeqMR(m, r, :) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_alpha_tmr));
                D_alpha_AeqMR(m, r, :) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_beta_tmr));
                D_beta_AeqMR(m, r, :) = sum(D_AeqMat,1);

                D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Mae0'*D_gamma_tmr));
                D_gamma_AeqMR(m, r, :) = sum(D_AeqMat,1);
                
                % BRM channel
                for b = 1:c.NB
                    trb = (c.B(:, b)-c.R(:, r))/norm(c.B(:, b)-c.R(:, r));
                    Psirb = c.Rae0'*c.RotR'*trb;            
                    tBFtemp = get_dir_from_angle(c.BFmatRB(:, r));  % get the BF direction vector for the b-th SA.
                    OmegaPsirb = Psirb - (c.Rae0'*tBFtemp)*c.beamsplit_coe;

                    trm = (c.M(:, m)-c.R(:, r))/norm(c.M(:, m)-c.R(:, r));
                    Psirm = c.Rae0'*c.RotR'*trm;            
                    tBFtemp = get_dir_from_angle(c.BFmatRM(:, r));  % get the BF direction vector for the b-th SA.
                    OmegaPsirm = Psirm - (c.Rae0'*tBFtemp)*c.beamsplit_coe;

                    AeqMat = exp(1j*2*pi./c.lambdak.*(OmegaPsirb + OmegaPsirm));
%                     AeqBRM(b, r, m, :) = real(sum(AeqMat,1));  % imaginary part is 0 after summation
                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*D_phiRM_trm));
                    D_phiRM_AeqBRM(b, r, m, :) = sum(D_AeqMat,1);

                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*D_thetaRM_trm));
                    D_thetaRM_AeqBRM(b, r, m, :) = sum(D_AeqMat,1);

                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*D_tauRM_trm));
                    D_tauRM_AeqBRM(b, r, m, :) = sum(D_AeqMat,1);
                   
                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*D_alpha_trm));
                    D_alpha_AeqBRM(b, r, m, :) = sum(D_AeqMat,1);

                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*D_beta_trm));
                    D_beta_AeqBRM(b, r, m, :) = sum(D_AeqMat,1);

                    D_AeqMat = AeqMat.*(1j*2*pi./c.lambdak.*(c.Rae0'*D_gamma_trm));
                    D_gamma_AeqBRM(b, r, m, :) = sum(D_AeqMat,1);
                end
            end
        end


    %     c.AeqMR = AeqMR;
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
    end

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