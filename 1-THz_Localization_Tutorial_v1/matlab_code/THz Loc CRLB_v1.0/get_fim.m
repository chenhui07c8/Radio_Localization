function [FIM0] = get_fim(c)
    % FIM matrix
    % 10x10 matrix for 2D
%     D_ub = zeros(1,10);
%     D_ub(1:4) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_tauBM_ub];
%     D_ub(5:8) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_tauRM_ub];
%     D_ub(9) = D_PhiM_ub(1);
%     D_ub(10) = D_beta_ub;

    % 14x14 matrix for 3D
%     D_ub = zeros(1,10);
%     D_ub(1:5) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_thetaBM_ub ,D_tauBM_ub];
%     D_ub(6:10) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_thetaRM_ub, D_tauRM_ub];
%     D_ub(11:13) = D_PhiM_ub;
%     D_ub(14) = D_beta_ub;

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
        
    c.D_PhiM_pm = get_rotation_over_position(c);
    
    % size(J)
    D_PhiM_pm = c.D_PhiM_pm;
    
    K = c.K;
    NB = c.NB;
    NM = c.NM;
    NR = c.NR;
    B = c.B;
    M = c.M;
    R = c.R;
    P = c.P;
    Xmk = c.Xmk;
    fk = c.fk;
    thetaBM = c.thetaBM;
    phiBM = c.phiBM;
    dBM = c.dBM;
    PM = c.PM;
    PB = c.PB;
    tBM = c.tBM;
    rhoBM = c.rhoBM;
    Omega = c.Omega;
    thetaRM = c.thetaRM;
    phiRM = c.phiRM;
    dRM = c.dRM;
    PR = c.PR;
    tRM = c.tRM;
    rhoBRM = c.rhoBRM;
    I_ns_nb = cell(K,NB);

    % AOSA parameters
    AeqBM = c.AeqBM;
    D_phiBM_AeqBM = c.D_phiBM_AeqBM;
    D_thetaBM_AeqBM = c.D_thetaBM_AeqBM;
    D_tauBM_AeqBM = c.D_tauBM_AeqBM;
    D_alpha_AeqBM = c.D_alpha_AeqBM;
    D_beta_AeqBM = c.D_beta_AeqBM;
    D_gamma_AeqBM = c.D_gamma_AeqBM;
    
    AeqBR = c.AeqBR;
    
    AeqMB = c.AeqMB;
    D_phiBM_AeqMB = c.D_phiBM_AeqMB;
    D_thetaBM_AeqMB = c.D_thetaBM_AeqMB;
    D_tauBM_AeqMB = c.D_tauBM_AeqMB;
    D_alpha_AeqMB = c.D_alpha_AeqMB;
    D_beta_AeqMB = c.D_beta_AeqMB;
    D_gamma_AeqMB = c.D_gamma_AeqMB;
    

    AeqMR = c.AeqMR;
    D_phiRM_AeqMR = c.D_phiRM_AeqMR;
    D_thetaRM_AeqMR = c.D_thetaRM_AeqMR;
    D_tauRM_AeqMR = c.D_tauRM_AeqMR;
    D_alpha_AeqMR = c.D_alpha_AeqMR;
    D_beta_AeqMR = c.D_beta_AeqMR;
    D_gamma_AeqMR = c.D_gamma_AeqMR;
    
    for k = 1:K
        % loop BS antennas
        for b = 1:NB
            % taubm = vecnorm(B(:, b)-M,2,1);
            ubm = zeros(NM, 1);
            D_thetaBM_taubm = zeros(NM, 1);
            D_phiBM_taubm = zeros(NM, 1);
            D_tauBM_taubm = zeros(NM, 1);
            ubrm = zeros(NM, NR);   % signal from the m-th MD antenna through r-th RIS to b-th BS antenna
            D_thetaRM_taurm = zeros(NM, NR);
            D_phiRM_taurm = zeros(NM, NR);
            D_tauRM_taurm = zeros(NM, NR);

            % Phi
            part1 = 0;
            D_PhiM_ubBM = zeros(NM, 3);
            part2 = 0;
            D_PhiM_ubBRM = zeros(NM, NR, 3);

            for m = 1:NM
                % BM *******
                taubm = norm(B(:,b) - M(:,m))/c.lightspeed;
                dbm = norm(B(:,b)-M(:,m));
                delta_bm = taubm-c.dBM/c.lightspeed;
%               ubm(m) = c.GantBM*Xmn(k, m)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*exp(-1j*2*pi*(c.fk(k)*taubm));
                ubm(m) = c.GantBM*Xmk(m, k)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*exp(-1j*2*pi*(c.fdk(k)*taubm + c.fc*delta_bm));
    %             dbm = taubmMat(b, m)*c;
    %             ubm(m) = Xmn(n, m)*exp(-1j*2*pi*fk(n)*(taubmMat(b, m)));
                % theta
                t_temp = [sind(thetaBM)*cosd(phiBM); sind(thetaBM)*sind(phiBM); -cosd(thetaBM)];
                D_thetaBM_taubm(m) = dBM/c.lightspeed/dbm*((M(:,m)-PM-B(:,b)+PB)'*t_temp);
                % phi
                t_temp = [cosd(thetaBM)*sind(phiBM); -cosd(thetaBM)*cosd(phiBM); 0];
                D_phiBM_taubm(m) = dBM/c.lightspeed/dbm*((M(:,m)-PM-B(:,b)+PB)'*t_temp);
                % tau
                Gbm2 = (M(:,m)-PM - B(:,b)+PB)'*tBM;
                D_tauBM_taubm(m) = 1/dbm*(dBM - Gbm2);

                % PhiM part 1
                D_PhiM_taubm = 1/c.lightspeed/dbm*((M(:,m)-PM-B(:,b)+PB)'*D_PhiM_pm{m} - dBM*(tBM'*D_PhiM_pm{m}));
                part1 = part1 + sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdac*c.xiBM)*ubm(m)*D_PhiM_taubm;
                D_PhiM_ubBM(m, :) = D_PhiM_taubm;
                % BRM ********
                for r = 1:NR
                    taubr = norm(B(:,b) - R(:,r))/c.lightspeed;
                    taurm = norm(R(:,r) - M(:,m))/c.lightspeed;
                    drm = norm(R(:,r)-M(:,m));
                    delta_brm = (taubr + taurm - c.dBR/c.lightspeed - c.dRM/c.lightspeed);
                    ubrm(m,r) = c.GantBRM*Xmk(m, k)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*Omega(r)*exp(-1j*2*pi*(c.fdk(k)*(taubr + taurm) + c.fc*delta_brm));
    %                 drm = taurmMat(r, m)*c;
    %                 ubrm(m,r) = Xmn(n, m)*Omega(n, r)*exp(-1j*2*pi*fk(n)*(taubrMat(b, r) + taurmMat(r, m)));
                    % theta
                    t_temp = [sind(thetaRM)*cosd(phiRM); sind(thetaRM)*sind(phiRM); -cosd(thetaRM)];
                    D_thetaRM_taurm(m,r) = dRM/c.lightspeed/drm*((M(:,m)-PM-R(:,r)+PR)'*t_temp);
                    % phi
                    t_temp = [cosd(thetaRM)*sind(phiRM); -cosd(thetaRM)*cosd(phiRM); 0];
                    D_phiRM_taurm(m,r) = dRM/c.lightspeed/drm*((M(:,m)-PM-R(:,r)+PR)'*t_temp);
                    % tau
                    Grm2 = (M(:,m)-PM - R(:,r)+PR)'*tRM;
                    D_tauRM_taurm(m,r) = 1/drm*(dRM - Grm2);

                    % PhiM part2
                    D_PhiM_taurm = 1/c.lightspeed/drm*((M(:,m)-PM-R(:,r)+PR)'*D_PhiM_pm{m} - dRM*(tRM'*D_PhiM_pm{m}));
                    part2 = part2 + sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdac*c.xiBM)*ubrm(m,r)*D_PhiM_taurm;
                    D_PhiM_ubBRM(m, r, :) = D_PhiM_taurm;
                end
            end

            
        
            
            % two-stage localization: added xi.
%             ubBM = sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdak(k)*c.xiBM)*sum(ubm, 1);
%             ubBRM = sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdak(k)*c.xiBRM)*sum(ubrm, [1, 2]);
% %             ubBM = sqrt(P)*rhoBM*sum(ubm);
% %             ubBRM = sqrt(P)*rhoBRM*sum(ubrm(:));
%             % BM
%             D_rhoBM_ub = ubBM/rhoBM;
%             D_xiBM_ub = 1j*2*pi/c.lambdak(k)*ubBM;
%             D_phiBM_ub = - 1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdak(k)*c.xiBM)*sum(ubm.*D_phiBM_taubm);
%             D_thetaBM_ub = -1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdak(k)*c.xiBM)*sum(ubm.*D_thetaBM_taubm);
%             D_tauBM_ub = - 1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdak(k)*c.xiBM)*sum(ubm.*D_tauBM_taubm);
% 
%             % BRM
%             D_rhoBRM_ub = ubBRM/rhoBRM;
%             D_xiBRM_ub = 1j*2*pi/c.lambdak(k)*ubBRM;
%             D_phiRM_ub   = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdak(k)*c.xiBRM)*sum(ubrm.*D_phiRM_taurm, 'all');
%             D_thetaRM_ub = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdak(k)*c.xiBRM)*sum(ubrm.*D_thetaRM_taurm, 'all');
%             D_tauRM_ub   = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdak(k)*c.xiBRM)*sum(ubrm.*D_tauRM_taurm, 'all');
% 
%             D_beta_ub = -1j*2*pi/c.lambdak(k)*(ubBM + ubBRM);
% 
%             D_PhiM_ub_1 = -1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdak(k)*c.xiBM)*sum(ubm.*D_PhiM_ubBM,1);
%             D_PhiM_ub_2 = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdak(k)*c.xiBRM)*reshape(sum(repmat(ubrm, [1 1 3]).*D_PhiM_ubBRM, [1 2]), [1,3]);
%             D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
%             D_PhiM_ub = -1j*2*pi*fk(k)*(part1 + part2);
            
% AOSA structure
            % sum over NM transmitter
            ubBM = sqrt(c.P)*rhoBM*exp(1j*2*pi/c.lambdac*c.xiBM)*sum(ubm.*AeqBM(b, :, k)'.*AeqMB(:, b, k), 1);
            ubBRM = sqrt(c.P)*rhoBRM*exp(1j*2*pi/c.lambdac*c.xiBRM)*sum(ubrm.*repmat(c.AeqBR(b, :, k), c.NM, 1).*AeqMR(:, :, k), [1, 2]);
% % 
            ubm = sqrt(P)*rhoBM*exp(1j*2*pi/c.lambdac*c.xiBM)*ubm;
            % BM channel (NM x 1)
            temp = AeqBM(b, :, k)'*-1j*2*pi*fk(k).*ubm.*D_phiBM_taubm.*AeqMB(:, b, k) + ...
                D_phiBM_AeqBM(b, :, k)'.*ubm.*AeqMB(:, b, k) + AeqBM(b, :, k)'.*ubm.*D_phiBM_AeqMB(:, b, k);
            D_phiBM_ub = sum(temp, 1);
            temp = AeqBM(b, :, k)'*-1j*2*pi*fk(k).*ubm.*D_thetaBM_taubm.*AeqMB(:, b, k) + ...
                D_thetaBM_AeqBM(b, :, k)'.*ubm.*AeqMB(:, b, k) + AeqBM(b, :, k)'.*ubm.*D_thetaBM_AeqMB(:, b, k);
            D_thetaBM_ub = sum(temp, 1);
            temp = AeqBM(b, :, k)'*-1j*2*pi*fk(k).*ubm.*D_tauBM_taubm.*AeqMB(:, b, k) + ...
                D_tauBM_AeqBM(b, :, k)'.*ubm.*AeqMB(:, b, k) + AeqBM(b, :, k)'.*ubm.*D_tauBM_AeqMB(:, b, k);
            D_tauBM_ub = sum(temp, 1);
%             D_tauBM_ub = -1j*2*pi*fk(k)*sum(AeqBM(b, :, k)'.*ubm.*D_tauBM_taubm.*AeqMB(:, b, k), 1);

            % BRM channel
            ubrm = sqrt(P)*rhoBRM*exp(1j*2*pi/c.lambdac*c.xiBRM)*ubrm;
            temp = AeqBR(b, :, k).*(-1j*2*pi*fk(k)*ubrm.*D_phiRM_taurm.*AeqMR(:, :, k) + ubrm.*D_phiRM_AeqMR(:, :, k));
            D_phiRM_ub   = sum(temp, [1, 2]);
            temp = AeqBR(b, :, k).*(-1j*2*pi*fk(k)*ubrm.*D_thetaRM_taurm.*AeqMR(:, :, k) + ubrm.*D_thetaRM_AeqMR(:, :, k));
            D_thetaRM_ub = sum(temp, [1, 2]);
            temp = AeqBR(b, :, k).*(-1j*2*pi*fk(k)*ubrm.*D_tauRM_taurm.*AeqMR(:, :, k) + ubrm.*D_tauRM_AeqMR(:, :, k));
            D_tauRM_ub = sum(temp, [1, 2]);
%             D_tauRM_ub   = -1j*2*pi*fk(k)*sum(AeqBR(b, :, k).*ubrm.*D_tauRM_taurm.*AeqMR(:, :,k), [1, 2]);

            % rotation
            D_PhiM_AeqBM = [D_alpha_AeqBM(b, :, k)', D_beta_AeqBM(b, :, k)', D_gamma_AeqBM(b, :, k)'];
            D_PhiM_AeqMB = [D_alpha_AeqMB(:, b, k), D_beta_AeqMB(:, b, k), D_gamma_AeqMB(:, b, k)];
            
            temp = AeqBM(b, :, k)'.*(-1j*2*pi*fk(k)*ubm.*D_PhiM_ubBM.*AeqMB(:, b, k) + ubm.*D_PhiM_AeqMB) + ...
                D_PhiM_AeqBM.*ubm.*AeqMB(:, b, k);
            D_PhiM_ub_1 = sum(temp, 1);
% 
            D_PhiM_AeqMR = zeros(NM, NR, 3);
            D_PhiM_AeqMR(:,:,1) = D_alpha_AeqMR(:, :, k);
            D_PhiM_AeqMR(:,:,2) = D_beta_AeqMR(:, :, k);
            D_PhiM_AeqMR(:,:,3) = D_gamma_AeqMR(:, :, k);
            
            temp = AeqBR(b, :, k).*(-1j*2*pi*fk(k)*repmat(ubrm, [1 1 3]).*D_PhiM_ubBRM.*AeqMR(:, :, k) + repmat(ubrm, [1 1 3]).*D_PhiM_AeqMR);
            D_PhiM_ub_2 = reshape(sum(temp, [1 2]), [1,3]);
          
            D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
            
            
%             D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
%             
            D_rhoBM_ub = ubBM/rhoBM;
            D_xiBM_ub = 1j*2*pi/c.lambdac*ubBM;
            D_rhoBRM_ub = ubBRM/rhoBRM;
            D_xiBRM_ub = 1j*2*pi/c.lambdac*ubBRM;
            D_beta_ub = -1j*2*pi/c.lambdak(k)*(ubBM + ubBRM);


            
%            3D done with calculations 
%             D_ub = zeros(1,14);
%             D_ub(1:5) = [D_rhoBM_ub,  D_xiBM_ub,  D_thetaBM_ub, D_phiBM_ub, D_tauBM_ub];
%             D_ub(6:10) = [D_rhoBRM_ub,D_xiBRM_ub, D_thetaRM_ub, D_phiRM_ub, D_tauRM_ub];
%             D_ub(11:13) = D_PhiM_ub;
%             D_ub(14) = D_beta_ub;
            
%             D_ub = zeros(1,10);
%             D_ub(1:4) = [D_phiBM_ub, D_tauBM_ub, D_phiRM_ub,  D_tauRM_ub]; % 2 3 5 6'
%             D_ub(5) = D_PhiM_ub(1);

%             D_ub(6:9) = [D_rhoBM_ub, D_xiBM_ub,  D_rhoBRM_ub, D_xiBRM_ub]; % 7-11 + 12
%             D_ub(10) = D_beta_ub;
%             D_ub(10) = D_PhiM_ub(1);
%             D_ub(1:4) = [D_phiBM_ub, D_tauBM_ub, D_phiRM_ub,  D_tauRM_ub]; % 2 3 5 6
%             D_ub(5:8) = [D_rhoBM_ub, D_xiBM_ub,  D_rhoBRM_ub, D_xiBRM_ub]; % 7-11 + 12
%             D_ub(9) = D_beta_ub;
%             D_ub(10) = D_PhiM_ub(1);
            if  c.pos_type == "2D"
            % 2D done with calculations 
                D_ub = zeros(1,10);
                D_ub(1:4) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_tauBM_ub];
                D_ub(5:8) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_tauRM_ub];
    %             D_ub(1:4) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_tauBM_ub/c.lightspeed];
    %             D_ub(5:8) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_tauRM_ub/c.lightspeed];
                D_ub(9) = D_PhiM_ub(1);
                D_ub(10) = D_beta_ub;
                if c.LR == 0
                    D_ub = D_ub([1:4 9:10]);
                end
            elseif  c.pos_type == "3D"
            % 3D done with calculations 
                D_ub = zeros(1,14);
                D_ub(1:5) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_thetaBM_ub ,D_tauBM_ub];
                D_ub(6:10) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_thetaRM_ub, D_tauRM_ub];
                D_ub(11:13) = D_PhiM_ub;
                D_ub(14) = D_beta_ub;
                if c.LR == 0
                    D_ub = D_ub([1:5 11:14]);
                end
            end
            I_ns_nb{k,b} = (D_ub'*D_ub);
        end
    end

    FIM0 = zeros(length(D_ub), length(D_ub));
    for k = 1:K
        I_ns_n = zeros(length(D_ub), length(D_ub));
        for b = 1:NB
            temp = I_ns_nb{k,b};
            I_ns_n = I_ns_n + temp;
        end
%         FIM = FIM + J*FIM0*J.';
        FIM0 = FIM0 + 2/c.sigma^2*real(I_ns_n);
    end
    
%     FIM0 = J*FIM0*J.';
    % 2D less unknowns
%     FIM = zeros(size(J,1), size(J,1));
%     for k = 1:K
%         I_ns_n = zeros(size(J,2), size(J,2));
%         for b = 1:NB
%             temp = I_ns_nb{k,b};
%             I_ns_n = I_ns_n + temp;
%         end
%         FIM0 = 2/c.sigma^2*real(I_ns_n);
%         FIM = FIM + J*FIM0*J.';
%     end

    % get crlb
%     if c.syn_type == "Syn"
%         CRLB = FIM(1:end-1,1:end-1)^(-1);
%     else
%         CRLB = FIM^(-1);
%     end

%     sqrt(diag(CRLB(1:3, 1:3)))

    
    
end