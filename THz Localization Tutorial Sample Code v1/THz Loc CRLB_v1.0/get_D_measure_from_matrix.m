

function c = get_D_measure_from_matrix(c)

    c.D_PhiM_pm = get_rotation_over_position(c);
    D_PhiM_pm = c.D_PhiM_pm;
    
    B = c.B;
    M = c.M;
    K = c.K;
    NB = c.NB;
    NM = c.NM;
    

    thetaBM = c.thetaBM;
    phiBM = c.phiBM;
    dBM = c.dBM;
    PM = c.PM;
    PB = c.PB;
    tBM = c.tBM;
%     rhoBM = c.rhoBM;


    % derivative matrix of BM channel
    D_thetaBM_taubm = zeros(NB,NM,K);
    D_phiBM_taubm = zeros(NB,NM,K);
    D_tauBM_taubm = zeros(NB,NM,K);
    D_PhiM_ubBM = zeros(NB,NM,K,3);
    
    
% LOS channel
    % calculate matrix derivative
    for k = 1:K
        for b = 1:NB
            for m = 1:NM
                % BM *******
                dbm = norm(B(:,b)-M(:,m));

                % theta
                t_temp = [sind(thetaBM)*cosd(phiBM); sind(thetaBM)*sind(phiBM); -cosd(thetaBM)];
                D_thetaBM_taubm(b, m, k) = dBM/c.lightspeed/dbm*((M(:,m)-PM-B(:,b)+PB)'*t_temp);
                % phi
                t_temp = [cosd(thetaBM)*sind(phiBM); -cosd(thetaBM)*cosd(phiBM); 0];
                D_phiBM_taubm(b, m, k) = dBM/c.lightspeed/dbm*((M(:,m)-PM-B(:,b)+PB)'*t_temp);
                % tau
                Gbm2 = (M(:,m)-PM - B(:,b)+PB)'*tBM;
                D_tauBM_taubm(b,m,k) = 1/dbm*(dBM - Gbm2);

                % PhiM part 1
                D_PhiM_ubBM(b,m,k,:) = 1/c.lightspeed/dbm*((M(:,m)-PM-B(:,b)+PB)'*D_PhiM_pm{m} - dBM*(tBM'*D_PhiM_pm{m}));
            end
        end
    end
    c.D_thetaBM_taubm = D_thetaBM_taubm;
    c.D_phiBM_taubm = D_phiBM_taubm;
    c.D_tauBM_taubm = D_tauBM_taubm;
    c.D_PhiM_ubBM = D_PhiM_ubBM;

    
% RIS channel 
    if (c.LR > 0)
        R = c.R;
        NR = c.NR;
        thetaRM = c.thetaRM;
        phiRM = c.phiRM;
        dRM = c.dRM;
        PR = c.PR;
        tRM = c.tRM;
        %     rhoBRM = c.rhoBRM;
        % derivative matrix of BRM channel
        D_thetaRM_taurm = zeros(NR,NM,K);
        D_phiRM_taurm = zeros(NR,NM,K);
        D_tauRM_taurm = zeros(NR,NM,K);
        D_PhiM_ubBRM = zeros(NR,NM,K,3);
        for k = 1:K
            for b = 1:NB
                for m = 1:NM
                    % BRM ********
                    for r = 1:NR
                        drm = norm(R(:,r)-M(:,m));
                        % theta
                        t_temp = [sind(thetaRM)*cosd(phiRM); sind(thetaRM)*sind(phiRM); -cosd(thetaRM)];
                        D_thetaRM_taurm(r,m,k) = dRM/c.lightspeed/drm*((M(:,m)-PM-R(:,r)+PR)'*t_temp);
                        % phi
                        t_temp = [cosd(thetaRM)*sind(phiRM); -cosd(thetaRM)*cosd(phiRM); 0];
                        D_phiRM_taurm(r,m,k) = dRM/c.lightspeed/drm*((M(:,m)-PM-R(:,r)+PR)'*t_temp);
                        % tau
                        Grm2 = (M(:,m)-PM - R(:,r)+PR)'*tRM;
                        D_tauRM_taurm(r,m,k) = 1/drm*(dRM - Grm2);

                        % PhiM part2
                        D_PhiM_ubBRM(r, m, k, :) = 1/c.lightspeed/drm*((M(:,m)-PM-R(:,r)+PR)'*D_PhiM_pm{m} - dRM*(tRM'*D_PhiM_pm{m}));
                    end
                end
            end
        end

        c.D_thetaRM_taurm = D_thetaRM_taurm;
        c.D_phiRM_taurm = D_phiRM_taurm;
        c.D_tauRM_taurm = D_tauRM_taurm;
        c.D_PhiM_ubBRM = D_PhiM_ubBRM;

    end
    
% NLOS channel
    if (c.LC>0)
        % derivative matrix of BC channel
        D_thetaBC_taubcm = zeros(NB,K,c.LC);
        D_phiBC_taubcm = zeros(NB,K,c.LC);
        D_tauBC_taubcm = zeros(NB,K,c.LC);
        D_thetaCM_taubcm = zeros(NM,K,c.LC);
        D_phiCM_taubcm = zeros(NM,K,c.LC);
        D_tauCM_taubcm = zeros(NM,K,c.LC);
        D_PhiM_ubBCM = zeros(NM,K,c.LC,3);  % only dependent on M
    
        % calculate matrix derivative
        for lc = 1:c.LC
            for k = 1:K
                for b = 1:NB
                    for m = 1:NM
                        % BC *******
                        dbc = norm(B(:,b)-c.PC(:,lc));
                        dmc = norm(c.PC(:,lc)-M(:,m));

                        % thetaBC
                        t_temp = [sind(c.thetaBC(lc))*cosd(c.phiBC(lc)); sind(c.thetaBC(lc))*sind(c.phiBC(lc)); -cosd(c.thetaBC(lc))];
                        D_thetaBC_taubcm(b, k, lc) = c.dBC(lc)/c.lightspeed/dbc*((-B(:,b)+PB)'*t_temp);
                        % phiBC
                        t_temp = [cosd(c.thetaBC(lc))*sind(c.phiBC(lc)); -cosd(c.thetaBC(lc))*cosd(c.phiBC(lc)); 0];
                        D_phiBC_taubcm(b, k, lc) = c.dBC(lc)/c.lightspeed/dbc*((-B(:,b)+PB)'*t_temp);
                        % tauBC
                        Gbc2 = (- B(:,b)+PB)'*c.tBC(:,lc);
                        D_tauBC_taubcm(b,k,lc) = 1/dbc*(c.dBC(lc) - Gbc2);

                        % thetaCM
                        t_temp = [sind(c.thetaCM(lc))*cosd(c.phiCM(lc)); sind(c.thetaCM(lc))*sind(c.phiCM(lc)); -cosd(c.thetaCM(lc))];
                        D_thetaCM_taubcm(m,k,lc) = c.dMC(lc)/c.lightspeed/dmc*((M(:,m)-PM)'*t_temp);
                        % phiCM
                        t_temp = [cosd(c.thetaCM(lc))*sind(c.phiCM(lc)); -cosd(c.thetaCM(lc))*cosd(c.phiCM(lc)); 0];
                        D_phiCM_taubcm(m,k,lc) = c.dMC(lc)/c.lightspeed/dmc*((M(:,m)-PM)'*t_temp);
                        % tauCM
                        Gcm2 = (M(:,m)-PM)'*(c.tCM(:,lc));
                        D_tauCM_taubcm(m,k,lc) = 1/dmc*(c.dMC(lc) - Gcm2);

                        % PhiM part3
                        D_PhiM_ubBCM(m, k, lc, :) = 1/c.lightspeed/dmc*((M(:,m)-PM)'*D_PhiM_pm{m} - c.dMC(lc)*(c.tCM(:, lc)'*D_PhiM_pm{m}));
                    end
                end
            end
        end
        c.D_thetaBC_taubcm = D_thetaBC_taubcm;
        c.D_phiBC_taubcm = D_phiBC_taubcm;
        c.D_tauBC_taubcm = D_tauBC_taubcm;
        c.D_thetaCM_taubcm = D_thetaCM_taubcm;
        c.D_phiCM_taubcm = D_phiCM_taubcm;
        c.D_tauCM_taubcm = D_tauCM_taubcm;
        c.D_PhiM_ubBCM = D_PhiM_ubBCM;  % only dependent on M
    
    end
    
end
