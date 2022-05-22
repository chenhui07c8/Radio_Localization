function [FIM_uplink, FIM_downlink] = get_fim_uplink(c)
    % FIM matrix
    % Jacobian matrix
%     c.J = get_jacobian_matrix(c);
%     J = c.J;

%     J = c.J([1 2 4], [1 3 4 5 7 8 9]);  % ignore theta, z, beta, gamma
%     J = c.J([1 2 4 7 8 9 10 11], [2 3 5 6 7 8 9 10 11 12]);  % ignore theta, z, beta, gamma
%     11 variables: x, y, z, a, b, gamma, phoBM, xiBM, rhoRM, xiRM, syn.
    % size(c.J)
    % Rotation over element position
    
    if c.AoSA_type == "SPP"
        c = get_D_sv_Aeq_SPP(c);
    elseif c.AoSA_type == "SSP"
        c = get_D_sv_Aeq_SSP(c);
    end
    
    % get derivatives of matrix
    c = get_D_measure_from_matrix(c);
    
    K = c.K;
    NB = c.NB;
    NM = c.NM;
%     NR = c.NR;
    P = c.P;
    fk = c.fk;
    I_ns_k_uplink = cell(K);
    I_ns_k_downlink = cell(K);

    
    % BM LOS channel
    if c.LB > 0
        rhoBM = c.rhoBM;

        D_thetaBM_taubm = c.D_thetaBM_taubm;
        D_phiBM_taubm = c.D_phiBM_taubm;
        D_tauBM_taubm = c.D_tauBM_taubm;
        D_PhiM_ubBM = c.D_PhiM_ubBM;
        
        % AOSA
        AeqBM = c.AeqBM;
        D_phiBM_AeqBM = c.D_phiBM_AeqBM;
        D_thetaBM_AeqBM = c.D_thetaBM_AeqBM;
        D_tauBM_AeqBM = c.D_tauBM_AeqBM;
        D_alpha_AeqBM = c.D_alpha_AeqBM;
        D_beta_AeqBM = c.D_beta_AeqBM;
        D_gamma_AeqBM = c.D_gamma_AeqBM;

        AeqMB = c.AeqMB;
        D_phiBM_AeqMB = c.D_phiBM_AeqMB;
        D_thetaBM_AeqMB = c.D_thetaBM_AeqMB;
        D_tauBM_AeqMB = c.D_tauBM_AeqMB;
        D_alpha_AeqMB = c.D_alpha_AeqMB;
        D_beta_AeqMB = c.D_beta_AeqMB;
        D_gamma_AeqMB = c.D_gamma_AeqMB;
    end
    
    % BRM RIS channel
    if c.LR > 0
        rhoBRM = c.rhoBRM;

        D_thetaRM_taurm = c.D_thetaRM_taurm;
        D_phiRM_taurm = c.D_phiRM_taurm;
        D_tauRM_taurm = c.D_tauRM_taurm;
        D_PhiM_ubBRM = c.D_PhiM_ubBRM;
        
        % AOSA
        AeqBR = c.AeqBR;
        AeqMR = c.AeqMR;
        D_phiRM_AeqMR = c.D_phiRM_AeqMR;
        D_thetaRM_AeqMR = c.D_thetaRM_AeqMR;
        D_tauRM_AeqMR = c.D_tauRM_AeqMR;
        D_alpha_AeqMR = c.D_alpha_AeqMR;
        D_beta_AeqMR = c.D_beta_AeqMR;
        D_gamma_AeqMR = c.D_gamma_AeqMR;
        
        % RIS side
        AeqBRM = c.AeqBRM;
        D_phiRM_AeqBRM = c.D_phiRM_AeqBRM;
        D_thetaRM_AeqBRM = c.D_thetaRM_AeqBRM;
        D_tauRM_AeqBRM = c.D_tauRM_AeqBRM;
        D_alpha_AeqBRM = c.D_alpha_AeqBRM;
        D_beta_AeqBRM = c.D_beta_AeqBRM;
        D_gamma_AeqBRM = c.D_gamma_AeqBRM;
    end
    
    % BCM NLOS channel
    if c.LC>0
        D_thetaBC_taubcm = c.D_thetaBC_taubcm;
        D_phiBC_taubcm = c.D_phiBC_taubcm;
        D_tauBC_taubcm = c.D_tauBC_taubcm;
        D_thetaCM_taubcm = c.D_thetaCM_taubcm;
        D_phiCM_taubcm = c.D_phiCM_taubcm;
        D_tauCM_taubcm = c.D_tauCM_taubcm;
        D_PhiM_ubBCM = c.D_PhiM_ubBCM;  % only dependent on M
        
        % BC channel
        AeqBC = c.AeqBC;
        D_phiBC_AeqBC = c.D_phiBC_AeqBC;
        D_thetaBC_AeqBC = c.D_thetaBC_AeqBC;
        D_tauBC_AeqBC = c.D_tauBC_AeqBC; % zeros in SPP
%         D_alpha_AeqBC = c.D_alpha_AeqBC; % zeros in SPP
%         D_beta_AeqBC = c.D_beta_AeqBC;  % zeros in SPP
%         D_gamma_AeqBC = c.D_gamma_AeqBC; % zeros in SPP

        % MC channel
        AeqMC = c.AeqMC;
        % AeqMC = zeros(c.NM, c.K, c.LC);
        D_phiCM_AeqMC = c.D_phiCM_AeqMC;
        D_thetaCM_AeqMC = c.D_thetaCM_AeqMC;
        D_tauCM_AeqMC = c.D_tauCM_AeqMC; % zeros in SPP
        D_alpha_AeqMC = c.D_alpha_AeqMC;
        D_beta_AeqMC = c.D_beta_AeqMC;
        D_gamma_AeqMC = c.D_gamma_AeqMC;
    end
    
    
    % calculate FIM
    for k = 1:K
        D_beta_ub = zeros(c.NB, 3);
        D_PhiM_ub_mat = zeros(c.NB, 3, 3);
        
        D_beta_um = zeros(c.NM, 3);
        D_PhiM_um_mat = zeros(c.NM, 3, 3);
        % BM channel 
        if (c.LB > 0)
            % uplink
            ubmMat = c.XiL(:,:,k);
            ubBM = c.uBM(:, k);
            D_rhoBM_ub = ubBM/rhoBM;
            D_xiBM_ub = 1j*2*pi/c.lambdac*ubBM;
            D_beta_ub(:, 1) = -1j*2*pi/c.lambdak(k)*(ubBM);
            
            umBM = c.uMB(:, k);
            D_rhoBM_um = umBM/rhoBM;
            D_xiBM_um = 1j*2*pi/c.lambdac*umBM;
            D_beta_um(:, 1) = -1j*2*pi/c.lambdak(k)*(umBM);

            % AOSA structure
            % BM channel
            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubmMat.*D_phiBM_taubm(:,:,k).*transpose(AeqMB(:,:,k)) + ...
                D_phiBM_AeqBM(:,:,k).*ubmMat.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubmMat.*transpose(D_phiBM_AeqMB(:,:,k));
            D_phiBM_ub = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
            D_phiBM_um = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);
            
            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubmMat.*D_thetaBM_taubm(:,:,k).*transpose(AeqMB(:,:,k)) + ...
                    D_thetaBM_AeqBM(:,:,k).*ubmMat.*transpose(AeqMB(:,:,k)) + ...
                    AeqBM(:,:,k).*ubmMat.*transpose(D_thetaBM_AeqMB(:,:,k));
            D_thetaBM_ub = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
            D_thetaBM_um = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

            temp = AeqBM(:,:,k).*(-1j*2*pi*fk(k)*ubmMat.*D_tauBM_taubm(:,:,k) + 1j*2*pi*c.fc*ubmMat).*transpose(AeqMB(:,:,k)) + ...
                    D_tauBM_AeqBM(:,:,k).*ubmMat.*transpose(AeqMB(:,:,k)) + ...
                    AeqBM(:,:,k).*ubmMat.*transpose(D_tauBM_AeqMB(:,:,k));
            D_tauBM_ub = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
            D_tauBM_um = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

            % orientation ub1
            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubmMat.*reshape(D_PhiM_ubBM(:,:,k,1),[NB NM]).*transpose(AeqMB(:,:,k)) + ...
                    D_alpha_AeqBM(:,:,k).*ubmMat.*transpose(AeqMB(:,:,k)) + ...
                    AeqBM(:,:,k).*ubmMat.*transpose(D_alpha_AeqMB(:,:,k));
            D_alpha_ub_1 = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
            D_alpha_um_1 = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubmMat.*reshape(D_PhiM_ubBM(:,:,k,2),[NB NM]).*transpose(AeqMB(:,:,k)) + ...
                    D_beta_AeqBM(:,:,k).*ubmMat.*transpose(AeqMB(:,:,k)) + ...
                    AeqBM(:,:,k).*ubmMat.*transpose(D_beta_AeqMB(:,:,k));
            D_beta_ub_1 = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
            D_beta_um_1 = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubmMat.*reshape(D_PhiM_ubBM(:,:,k,3),[NB NM]).*transpose(AeqMB(:,:,k)) + ...
                    D_gamma_AeqBM(:,:,k).*ubmMat.*transpose(AeqMB(:,:,k)) + ...
                    AeqBM(:,:,k).*ubmMat.*transpose(D_gamma_AeqMB(:,:,k));
            D_gamma_ub_1 = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
            D_gamma_um_1 = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

            D_PhiM_ub_mat(:,:,1) = [D_alpha_ub_1 D_beta_ub_1 D_gamma_ub_1];
            D_PhiM_um_mat(:,:,1) = [D_alpha_um_1 D_beta_um_1 D_gamma_um_1];
        end
        

        % BRM channel
        if (c.LR > 0)
            
            ubrMat = c.XiBR(:,:,k);
            urmMat = c.XiRM(:,:,k);
            ubBRM = c.uBRM(:, k);
            D_rhoBRM_ub = ubBRM/rhoBRM;
            D_xiBRM_ub = 1j*2*pi/c.lambdac*ubBRM;
            D_beta_ub(:, 2) = -1j*2*pi/c.lambdak(k)*(ubBRM);
            
            umBRM = c.uMRB(:, k);
            D_rhoBRM_um = umBRM/rhoBRM;
            D_xiBRM_um = 1j*2*pi/c.lambdac*umBRM;
            D_beta_um(:, 2) = -1j*2*pi/c.lambdak(k)*(umBRM);

            D_phiRM_ubm = zeros(c.NB, c.NM);
            D_thetaRM_ubm = zeros(c.NB, c.NM);
            D_tauRM_ubm = zeros(c.NB, c.NM);
            D_alpha_ubm = zeros(c.NB, c.NM);
            D_beta_ubm = zeros(c.NB, c.NM);
            D_gamma_ubm = zeros(c.NB, c.NM);

            if(c.NsaR > 1)  % AOSA at RIS
                for b = 1:c.NB
                    for m = 1:c.NM
                        ubr = c.XiBR(b,:,k);
                        urm = c.XiRM(:,m,k);
                        % AeqBR(b,r,k)*(ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*urm)*AeqMR(m,r)
                        % urm' + AeqMR' + AeqBRM'
                        % phiRM

                        D_phiRM_ubm(b, m) = AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))...
                            *(urm*-1j*2*pi*fk(k).*D_phiRM_taurm(:,m,k).*transpose(AeqMR(m,:,k))) + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(D_phiRM_AeqBRM(b, :, m, k)))*(urm.*(AeqMR(m,:,k)).') + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(urm.*(D_phiRM_AeqMR(m,:,k)).');

                        D_thetaRM_ubm(b, m) = AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))...
                            *(urm*-1j*2*pi*fk(k).*D_thetaRM_taurm(:,m,k).*transpose(AeqMR(m,:,k))) + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(D_thetaRM_AeqBRM(b, :, m, k)))*(urm.*(AeqMR(m,:,k)).') + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(urm.*(D_thetaRM_AeqMR(m,:,k)).');

                        D_tauRM_ubm(b, m) = AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))...
                            *((urm*-1j*2*pi*fk(k).*D_tauRM_taurm(:,m,k)+1j*2*pi*c.fc*urm).*transpose(AeqMR(m,:,k))) + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(D_tauRM_AeqBRM(b, :, m, k)))*(urm.*(AeqMR(m,:,k)).') + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(urm.*(D_tauRM_AeqMR(m,:,k)).');

                        D_alpha_ubm(b, m) = AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))...
                            *(urm*-1j*2*pi*fk(k).*D_PhiM_ubBRM(: ,m ,k, 1).*transpose(AeqMR(m,:,k))) + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(D_alpha_AeqBRM(b, :, m, k)))*(urm.*(AeqMR(m,:,k)).') + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(urm.*(D_alpha_AeqMR(m,:,k)).');

                        D_beta_ubm(b, m) = AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))...
                            *(urm*-1j*2*pi*fk(k).*D_PhiM_ubBRM(:, m ,k, 2).*transpose(AeqMR(m,:,k))) + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(D_beta_AeqBRM(b, :, m, k)))*(urm.*(AeqMR(m,:,k)).') + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(urm.*(D_beta_AeqMR(m,:,k)).');

                        D_gamma_ubm(b, m) = AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))...
                            *(urm*-1j*2*pi*fk(k).*D_PhiM_ubBRM(:, m ,k, 3).*transpose(AeqMR(m,:,k))) + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(D_gamma_AeqBRM(b, :, m, k)))*(urm.*(AeqMR(m,:,k)).') + ...
                            AeqBR(b,:,k).*ubr*diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(urm.*(D_gamma_AeqMR(m,:,k)).');
                    end
                end
                % uplink
                D_phiRM_ub = sqrt(P)*c.alphaR(k)*D_phiRM_ubm*c.Xmk(:,k);
                D_thetaRM_ub = sqrt(P)*c.alphaR(k)*D_thetaRM_ubm*c.Xmk(:,k);
                D_tauRM_ub = sqrt(P)*c.alphaR(k)*D_tauRM_ubm*c.Xmk(:,k);

                D_alpha_ub_2 = sqrt(P)*c.alphaR(k)*D_alpha_ubm*c.Xmk(:,k);
                D_beta_ub_2 = sqrt(P)*c.alphaR(k)*D_beta_ubm*c.Xmk(:,k);
                D_gamma_ub_2 = sqrt(P)*c.alphaR(k)*D_gamma_ubm*c.Xmk(:,k);
                D_PhiM_ub_mat(:,:,2) = [D_alpha_ub_2 D_beta_ub_2 D_gamma_ub_2];
                
                % downlink
                D_phiRM_um = sqrt(P)*c.alphaR(k)*D_phiRM_ubm.'*c.Xbk(:,k);
                D_thetaRM_um = sqrt(P)*c.alphaR(k)*D_thetaRM_ubm.'*c.Xbk(:,k);
                D_tauRM_um = sqrt(P)*c.alphaR(k)*D_tauRM_ubm.'*c.Xbk(:,k);

                D_alpha_um_2 = sqrt(P)*c.alphaR(k)*D_alpha_ubm.'*c.Xbk(:,k);
                D_beta_um_2 = sqrt(P)*c.alphaR(k)*D_beta_ubm.'*c.Xbk(:,k);
                D_gamma_um_2 = sqrt(P)*c.alphaR(k)*D_gamma_ubm.'*c.Xbk(:,k);
                D_PhiM_um_mat(:,:,2) = [D_alpha_um_2 D_beta_um_2 D_gamma_um_2];
            else % no AOSA at RIS (should be faster..)
                temp = AeqBR(:,:,k).*ubrMat*diag(c.Omega)* ...
                    (-1j*2*pi*fk(k)*urmMat.*D_phiRM_taurm(:,:,k).*AeqMR(:,:,k).' + urmMat.*D_phiRM_AeqMR(:,:,k).');
                D_phiRM_ub = sqrt(P)*c.alphaR(k)*temp*c.Xmk(:,k);
                D_phiRM_um = sqrt(P)*c.alphaR(k)*temp.'*c.Xbk(:,k);

                temp = AeqBR(:,:,k).*ubrMat*diag(c.Omega)*...
                    (-1j*2*pi*fk(k)*urmMat.*D_thetaRM_taurm(:,:,k).*AeqMR(:, :, k).' + urmMat.*D_thetaRM_AeqMR(:,:,k).');
                D_thetaRM_ub = sqrt(P)*c.alphaR(k)*temp*c.Xmk(:,k);
                D_thetaRM_um = sqrt(P)*c.alphaR(k)*temp.'*c.Xbk(:,k);

                temp = AeqBR(:,:,k).*ubrMat*diag(c.Omega)*...
                    ((-1j*2*pi*fk(k)*urmMat.*D_tauRM_taurm(:,:,k) + 1j*2*pi*c.fc*urmMat).*AeqMR(:,:,k).' + urmMat.*D_tauRM_AeqMR(:,:,k).');
                D_tauRM_ub = sqrt(P)*c.alphaR(k)*temp*c.Xmk(:,k);
                D_tauRM_um = sqrt(P)*c.alphaR(k)*temp.'*c.Xbk(:,k);

%                 orientation ub2
                temp = AeqBR(:,:,k).*ubrMat*diag(c.Omega)*...
                    (-1j*2*pi*fk(k)*urmMat.*D_PhiM_ubBRM(:,:,k,1).*AeqMR(:,:,k).' + urmMat.*D_alpha_AeqMR(:,:,k).');
                D_alpha_ub_2 = sqrt(P)*c.alphaR(k)*temp*c.Xmk(:,k);
                D_alpha_um_2 = sqrt(P)*c.alphaR(k)*temp.'*c.Xbk(:,k);

                temp = AeqBR(:,:,k).*ubrMat*diag(c.Omega)*...
                    (-1j*2*pi*fk(k)*urmMat.*D_PhiM_ubBRM(:,:,k,2).*AeqMR(:,:,k).' + urmMat.*D_beta_AeqMR(:,:,k).');
                D_beta_ub_2 = sqrt(P)*c.alphaR(k)*temp*c.Xmk(:,k);
                D_beta_um_2 = sqrt(P)*c.alphaR(k)*temp.'*c.Xbk(:,k);

                temp = AeqBR(:,:,k).*ubrMat*diag(c.Omega)*...
                    (-1j*2*pi*fk(k)*urmMat.*D_PhiM_ubBRM(:,:,k,3).*AeqMR(:,:,k).' + urmMat.*D_gamma_AeqMR(:,:,k).');
                D_gamma_ub_2 = sqrt(P)*c.alphaR(k)*temp*c.Xmk(:,k);
                D_gamma_um_2 = sqrt(P)*c.alphaR(k)*temp.'*c.Xbk(:,k);

                D_PhiM_ub_mat(:,:,2) = [D_alpha_ub_2 D_beta_ub_2 D_gamma_ub_2];
                D_PhiM_um_mat(:,:,2) = [D_alpha_um_2 D_beta_um_2 D_gamma_um_2];
            end
        end
        
        if (c.LC>0)
            D_rhoBCM_ub = zeros(c.NB, c.LC);
            D_xiBCM_ub = zeros(c.NB, c.LC);
            
            D_phiBC_ub = zeros(c.NB, c.LC);
            D_thetaBC_ub = zeros(c.NB, c.LC);
            D_tauBC_ub = zeros(c.NB, c.LC);
            D_phiCM_ub = zeros(c.NB, c.LC);
            D_thetaCM_ub = zeros(c.NB, c.LC);
            D_tauCM_ub = zeros(c.NB, c.LC);
            
            D_rhoBCM_um = zeros(c.NM, c.LC);
            D_xiBCM_um = zeros(c.NM, c.LC);
            
            D_phiBC_um = zeros(c.NM, c.LC);
            D_thetaBC_um = zeros(c.NM, c.LC);
            D_tauBC_um = zeros(c.NM, c.LC);
            D_phiCM_um = zeros(c.NM, c.LC);
            D_thetaCM_um = zeros(c.NM, c.LC);
            D_tauCM_um = zeros(c.NM, c.LC);
%             D_PhiM_ub_3 = zeros(c.NB, c.LC, 3);
            
            % clock offset beta
            ubBCM = c.uBCM(:,k);
            D_beta_ub(:, 3) = D_beta_ub(:, 3) -1j*2*pi/c.lambdak(k)*ubBCM;
            
            umBCM = c.uMCB(:,k);
            D_beta_um(:, 3) = D_beta_um(:, 3) -1j*2*pi/c.lambdak(k)*umBCM;

            for l = 1:c.LC
                ubcMat = c.XiBC(:,k,l);
                ucmMat = c.XiMC(:,k,l);
                
                ubBCML = c.uBCML(:,k,l);
                D_rhoBCM_ub(:, l) = ubBCML/c.rhoBCM(l);
                D_xiBCM_ub(:, l) = 1j*2*pi/c.lambdac*ubBCML;
                umMCBL = c.uMCBL(:,k,l);
                D_rhoBCM_um(:, l) = umMCBL/c.rhoBCM(l);
                D_xiBCM_um(:, l) = 1j*2*pi/c.lambdac*umMCBL;
                
                % AOSA structure
                % BC channel
                temp = AeqBC(:,k,l)*-1j*2*pi*fk(k).*ubcMat.*D_phiBC_taubcm(:,k,l).*transpose(ucmMat.*AeqMC(:,k,l)) + ...
                    D_phiBC_AeqBC(:,k,l).*ubcMat.*transpose(ucmMat.*AeqMC(:,k,l));
                D_phiBC_ub(:,l) = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_phiBC_um(:,l) = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                temp = AeqBC(:,k,l)*-1j*2*pi*fk(k).*ubcMat.*D_thetaBC_taubcm(:,k,l).*transpose(ucmMat.*AeqMC(:,k,l)) + ...
                        D_thetaBC_AeqBC(:,k,l).*ubcMat.*transpose(ucmMat.*AeqMC(:,k,l));
                D_thetaBC_ub(:,l) = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_thetaBC_um(:,l) = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                temp = AeqBC(:,k,l).*(-1j*2*pi*fk(k)*ubcMat.*D_tauBC_taubcm(:,k,l) + 1j*2*pi*c.fc*ubcMat).*transpose(ucmMat.*AeqMC(:,k,l)) + ...
                        D_tauBC_AeqBC(:,k,l).*ubcMat.*transpose(ucmMat.*AeqMC(:,k,l));
                D_tauBC_ub(:,l) = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_tauBC_um(:,l) = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                
                % CM channel
                temp = AeqBC(:,k,l).*ubcMat* ...
                    transpose(-1j*2*pi*fk(k)*ucmMat.*D_phiCM_taubcm(:,k,l).*AeqMC(:,k,l) + ucmMat.*D_phiCM_AeqMC(:,k,l));
                D_phiCM_ub(:,l) = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_phiCM_um(:,l) = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                temp = AeqBC(:,k,l).*ubcMat*...
                    transpose(-1j*2*pi*fk(k)*ucmMat.*D_thetaCM_taubcm(:,k,l).*AeqMC(:,k,l) + ucmMat.*D_thetaCM_AeqMC(:,k,l));
                D_thetaCM_ub(:,l) = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_thetaCM_um(:,l) = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                temp = AeqBC(:,k,l).*ubcMat*...
                    transpose((-1j*2*pi*fk(k)*ucmMat.*D_tauCM_taubcm(:,k,l) + 1j*2*pi*c.fc*ucmMat).*AeqMC(:,k,l) + ucmMat.*D_tauCM_AeqMC(:,k,l));
                D_tauCM_ub(:,l) = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_tauCM_um(:,l) = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                % orientation ub3
                temp = AeqBC(:,k,l).*ubcMat*...
                    transpose(-1j*2*pi*fk(k)*ucmMat.*reshape(D_PhiM_ubBCM(:,k,l,1),[NM 1]).*AeqMC(:,k,l) + ucmMat.*D_alpha_AeqMC(:,k,l));
                D_alpha_ub_3 = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_alpha_um_3 = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                temp = AeqBC(:,k,l).*ubcMat*...
                    transpose(-1j*2*pi*fk(k)*ucmMat.*reshape(D_PhiM_ubBCM(:,k,l,2),[NM 1]).*AeqMC(:,k,l) + ucmMat.*D_beta_AeqMC(:,k,l));
                D_beta_ub_3 = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_beta_um_3 = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                temp = AeqBC(:,k,l).*ubcMat*...
                    transpose(-1j*2*pi*fk(k)*ucmMat.*reshape(D_PhiM_ubBCM(:,k,l,3),[NM 1]).*AeqMC(:,k,l) + ucmMat.*D_gamma_AeqMC(:,k,l));
                D_gamma_ub_3 = sqrt(P)*c.alphaN(k,l)*temp*c.Xmk(:,k);
                D_gamma_um_3 = sqrt(P)*c.alphaN(k,l)*temp.'*c.Xbk(:,k);

                D_PhiM_ub_mat(:,:,3) = D_PhiM_ub_mat(:,:,3) + [D_alpha_ub_3 D_beta_ub_3 D_gamma_ub_3];
                D_PhiM_um_mat(:,:,3) = D_PhiM_um_mat(:,:,3) + [D_alpha_um_3 D_beta_um_3 D_gamma_um_3];

            end
        end
        D_PhiM_ub = sum(D_PhiM_ub_mat, 3);
        D_beta_ub = sum(D_beta_ub, 2);
        
        D_PhiM_um = sum(D_PhiM_um_mat, 3);
        D_beta_um = sum(D_beta_um, 2);
        

% uplink %%%%%%%%%%%%%%%%%%%%%%%%%
        D_ub = zeros(c.NB, c.LB*5 + c.LR*5 + c.LC*8 + 3 + 1);
        if (c.LB > 0)
            for l = 1:c.LB
            	D_ub(:, (1:5) + 5*(l-1)) = [D_rhoBM_ub(:, l),  D_xiBM_ub(:, l),  D_phiBM_ub(:, l),...
                    D_thetaBM_ub(:, l), D_tauBM_ub(:, l)];
            end
        end
        if (c.LR > 0)
            for l = 1:c.LR
            	D_ub(:, 5*c.LB + (1:5) + 5*(l-1)) = [D_rhoBRM_ub(:, l),  D_xiBRM_ub(:, l),  D_phiRM_ub(:, l),...
                    D_thetaRM_ub(:, l), D_tauRM_ub(:, l)];
            end
        end
        D_ub(:, end-3:end-1) = D_PhiM_ub;
        D_ub(:, end) = D_beta_ub;
        if (c.LC > 0)
            % each cluster has 6 unknowns in 2D (rho, xi, phiBC, tauBC, phiCM, tauCM)
            for l = 1:c.LC
                D_ub(:, (5*c.LB + 5*c.LR +(1:8) + 8*(l-1))) = [D_rhoBCM_ub(:,l),  D_xiBCM_ub(:,l),  D_phiBC_ub(:,l), D_thetaBC_ub(:,l), ...
                    D_tauBC_ub(:,l), D_phiCM_ub(:,l), D_thetaCM_ub(:,l), D_tauCM_ub(:,l)];
            end
        end
        I_ns_k_uplink{k} = (D_ub'*D_ub); 
        
% downlink %%%%%%%%%%%%%%%%%%%%%%%%
        D_um = zeros(c.NM, c.LB*5 + c.LR*5 + c.LC*8 + 3 + 1);
        if (c.LB > 0)
            for l = 1:c.LB
            	D_um(:, (1:5) + 5*(l-1)) = [D_rhoBM_um(:, l),  D_xiBM_um(:, l),  D_phiBM_um(:, l),...
                    D_thetaBM_um(:, l), D_tauBM_um(:, l)];
            end
        end
        if (c.LR > 0)
            for l = 1:c.LR
            	D_um(:, 5*c.LB + (1:5) + 5*(l-1)) = [D_rhoBRM_um(:, l),  D_xiBRM_um(:, l),  D_phiRM_um(:, l),...
                    D_thetaRM_um(:, l), D_tauRM_um(:, l)];
            end
        end
        D_um(:, end-3:end-1) = D_PhiM_um;
        D_um(:, end) = D_beta_um;
        if (c.LC > 0)
            % each cluster has 6 unknowns in 2D (rho, xi, phiBC, tauBC, phiCM, tauCM)
            for l = 1:c.LC
                D_um(:, (5*c.LB + 5*c.LR +(1:8) + 8*(l-1))) = [D_rhoBCM_um(:,l),  D_xiBCM_um(:,l),  ...
                    D_phiBC_um(:,l), D_thetaBC_um(:,l), D_tauBC_um(:,l), ...
                    D_phiCM_um(:,l), D_thetaCM_um(:,l), D_tauCM_um(:,l)];
            end
        end
        I_ns_k_downlink{k} = (D_um'*D_um);
    end

    
    FIM0 = zeros(size(D_ub,2), size(D_ub,2));
    for k = 1:K
        I_ns_n = zeros(size(FIM0));
        temp = I_ns_k_uplink{k};
        I_ns_n = I_ns_n + temp;
        FIM0 = FIM0 + 2/c.sigma^2*real(I_ns_n);
    end
    FIM_uplink = FIM0;
    
    FIM0 = zeros(size(D_um,2), size(D_um,2));
    for k = 1:K
        I_ns_n = zeros(size(FIM0));
        temp = I_ns_k_downlink{k};
        I_ns_n = I_ns_n + temp;
        FIM0 = FIM0 + 2/c.sigma^2*real(I_ns_n);
    end
    FIM_downlink = FIM0;
    
    
end