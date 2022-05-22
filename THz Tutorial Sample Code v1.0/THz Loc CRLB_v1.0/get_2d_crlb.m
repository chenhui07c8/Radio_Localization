function [CRLB, FIM] = get_2d_crlb(c)
    % Jacobian matrix
    c.J = get_jacobian_matrix(c);
%     J = p.J([1 2 4], [1 3 4 5 7 8 9]);  % ignore theta, z, beta, gamma
%     J = p.J([1 2 4 7 8 9 10 11], [2 3 5 6 7 8 9 10 11 12]);  % ignore theta, z, beta, gamma
%     11 variables: x, y, z, a, b, gamma, phoBM, xiBM, rhoRM, xiRM, syn.
    % size(p.J)
    % Rotation over element position
    c.D_PhiM_pm = get_rotation_over_position(c);
    
    J = c.J;

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
%     c = p.c;
    Xmn = c.Xmn;
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
    D_phi_AeqBM = c.D_phi_AeqBM;
    D_theta_AeqBM = c.D_theta_AeqBM;
    
    AeqBR = c.AeqBR;
    
    AeqMB = c.AeqMB;
    D_phi_AeqMB = c.D_phi_AeqMB;
    D_theta_AeqMB = c.D_theta_AeqMB;
    D_alpha_AeqMB = c.D_alpha_AeqMB;
    D_beta_AeqMB = c.D_beta_AeqMB;
    D_gamma_AeqMB = c.D_gamma_AeqMB;
    

    AeqMR = c.AeqMR;
    D_phi_AeqMR = c.D_phi_AeqMR;
    D_theta_AeqMR = c.D_theta_AeqMR;
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
%               ubm(m) = p.GantBM*Xmn(k, m)*exp(-1j*2*pi/p.lambdak(k)*p.beta)*exp(-1j*2*pi*(p.fk(k)*taubm));
                ubm(m) = c.GantBM*Xmn(k, m)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*exp(-1j*2*pi*(c.fdk(k)*taubm + c.fc*delta_bm));
    %             dbm = taubmMat(b, m)*p.c;
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
                part1 = part1 + sqrt(P)*rhoBM*exp(1j*c.xiBM)*ubm(m)*D_PhiM_taubm;
                D_PhiM_ubBM(m, :) = D_PhiM_taubm;
                % BRM ********
                for r = 1:NR
                    taubr = norm(B(:,b) - R(:,r))/c.lightspeed;
                    taurm = norm(R(:,r) - M(:,m))/c.lightspeed;
                    drm = norm(R(:,r)-M(:,m));
                    delta_brm = (taubr + taurm - c.dBR/c.lightspeed - c.dRM/c.lightspeed);
                    ubrm(m,r) = c.GantBRM*Xmn(k, m)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*Omega(r)*exp(-1j*2*pi*(c.fdk(k)*(taubr + taurm) + c.fc*delta_brm));
    %                 drm = taurmMat(r, m)*p.c;
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
                    part2 = part2 + sqrt(P)*rhoBRM*exp(1j*c.xiBM)*ubrm(m,r)*D_PhiM_taurm;
                    D_PhiM_ubBRM(m, r, :) = D_PhiM_taurm;
                end
            end

            
        
            
            % two-stage localization: added xi.
%             ubBM = sqrt(P)*rhoBM*exp(1j*p.xiBM)*sum(ubm, 1);
%             ubBRM = sqrt(P)*rhoBRM*exp(1j*p.xiBRM)*sum(ubrm, [1, 2]);
% %             ubBM = sqrt(P)*rhoBM*sum(ubm);
% %             ubBRM = sqrt(P)*rhoBRM*sum(ubrm(:));
%             % BM
%             D_rhoBM_ub = ubBM/rhoBM;
%             D_xiBM_ub = 1j*ubBM;
%             D_phiBM_ub = - 1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*p.xiBM)*sum(ubm.*D_phiBM_taubm);
%             D_thetaBM_ub = -1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*p.xiBM)*sum(ubm.*D_thetaBM_taubm);
%             D_tauBM_ub = - 1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*p.xiBM)*sum(ubm.*D_tauBM_taubm);
% 
%             % BRM
%             D_rhoBRM_ub = ubBRM/rhoBRM;
%             D_xiBRM_ub = 1j*ubBRM;
%             D_phiRM_ub   = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*p.xiBRM)*sum(ubrm.*D_phiRM_taurm, 'all');
%             D_thetaRM_ub = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*p.xiBRM)*sum(ubrm.*D_thetaRM_taurm, 'all');
%             D_tauRM_ub   = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*p.xiBRM)*sum(ubrm.*D_tauRM_taurm, 'all');
% 
%             D_beta_ub = -1j*2*pi/p.lambdak(k)*(ubBM + ubBRM);
% 
%             D_PhiM_ub_1 = -1j*2*pi*fk(k)*sqrt(P)*rhoBM*exp(1j*p.xiBM)*sum(ubm.*D_PhiM_ubBM,1);
%             D_PhiM_ub_2 = -1j*2*pi*fk(k)*sqrt(P)*rhoBRM*exp(1j*p.xiBRM)*reshape(sum(repmat(ubrm, [1 1 3]).*D_PhiM_ubBRM, [1 2]), [1,3]);
%             D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
%             D_PhiM_ub = -1j*2*pi*fk(k)*(part1 + part2);
            
% AOSA structure
            % sum over NM transmitter
            ubBM = sqrt(c.P)*rhoBM*exp(1j*c.xiBM)*c.AeqBM(b,k)*sum(ubm.*AeqMB(:,k), 1);
            ubBRM = sqrt(c.P)*rhoBRM*exp(1j*c.xiBRM)*c.AeqBR(b,k)*sum(ubrm.*AeqMR(:,k), [1, 2]);
% % 
            ubm = sqrt(P)*rhoBM*exp(1j*c.xiBM)*ubm;
            % BM channel
            temp = AeqBM(b, k)*-1j*2*pi*fk(k)*ubm.*D_phiBM_taubm.*AeqMB(:, k) + D_phi_AeqBM(b, k)*ubm.*AeqMB(:,k) + AeqBM(b, k)*ubm.*D_phi_AeqMB(:,k);
            D_phiBM_ub = sum(temp, 1);
            temp = AeqBM(b, k)*-1j*2*pi*fk(k)*ubm.*D_thetaBM_taubm.*AeqMB(:, k) + D_theta_AeqBM(b, k)*ubm.*AeqMB(:,k) + AeqBM(b, k)*ubm.*D_theta_AeqMB(:,k);
            D_thetaBM_ub = sum(temp, 1);
            D_tauBM_ub = c.AeqBM(b,k)*-1j*2*pi*fk(k)*sum(ubm.*D_tauBM_taubm.*AeqMB(:,k), 1);

%             BRM channel
            ubrm = sqrt(P)*rhoBRM*exp(1j*c.xiBRM)*ubrm;
%           temp = AeqBR(b, k)*(-1j*2*pi*fk(k)*ubrm.*D_phiRM_taurm.*AeqMR(:,k) + ubrm.*D_phi_AeqMR(:, k))/sqrt(P)/rhoBRM/exp(1j*p.xiBRM);
            temp = AeqBR(b, k)*(-1j*2*pi*fk(k)*ubrm.*D_phiRM_taurm.*AeqMR(:,k) + ubrm.*D_phi_AeqMR(:, k));
            D_phiRM_ub   = sum(temp, [1, 2]);
            temp = AeqBR(b, k)*(-1j*2*pi*fk(k)*ubrm.*D_thetaRM_taurm.*AeqMR(:,k) + ubrm.*D_theta_AeqMR(:, k));
            D_thetaRM_ub = sum(temp, [1, 2]);
            D_tauRM_ub   = c.AeqBR(b,k)*-1j*2*pi*fk(k)*sum(ubrm.*D_tauRM_taurm.*AeqMR(:,k), [1, 2]);

            % rotation
            D_PhiM_AeqMB = [D_alpha_AeqMB(:, k), D_beta_AeqMB(:, k), D_gamma_AeqMB(:, k)];
            temp = AeqBM(b, k)*(-1j*2*pi*fk(k)*ubm.*D_PhiM_ubBM.*AeqMB(:, k) + ubm.*D_PhiM_AeqMB);
            D_PhiM_ub_1 = sum(temp, 1);
% 
            D_PhiM_AeqMR = zeros(NM, NR, 3);
            D_PhiM_AeqMR(:,:,1) = repmat(D_alpha_AeqMR(:, k), [1 NR]);
            D_PhiM_AeqMR(:,:,2) = repmat(D_beta_AeqMR(:, k), [1 NR]);
            D_PhiM_AeqMR(:,:,3) = repmat(D_gamma_AeqMR(:, k), [1 NR]);
            temp = AeqBR(b, k)*(-1j*2*pi*fk(k)*repmat(ubrm, [1 1 3]).*D_PhiM_ubBRM.*AeqMB(:, k) + repmat(ubrm, [1 1 3]).*D_PhiM_AeqMR);
            D_PhiM_ub_2 = reshape(sum(temp, [1 2]), [1,3]);
          
            D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
            
            
%             D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
%             
            D_rhoBM_ub = ubBM/rhoBM;
            D_xiBM_ub = 1j*ubBM;
            D_rhoBRM_ub = ubBRM/rhoBRM;
            D_xiBRM_ub = 1j*ubBRM;
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

            % 2D done with calculations 
%             D_ub = zeros(1,10);
            D_ub(1:4) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_tauBM_ub];
            D_ub(5:8) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_tauRM_ub];
%             D_ub(1:4) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_tauBM_ub/p.c];
%             D_ub(5:8) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_tauRM_ub/p.c];
            D_ub(9) = D_PhiM_ub(1);
            D_ub(10) = D_beta_ub;

            I_ns_nb{k,b} = (D_ub'*D_ub);
        end
    end

    % 2D less unknowns
    FIM = zeros(size(J,2), size(J,2));
    for k = 1:K
        I_ns_n = zeros(size(J,2), size(J,2));
        for b = 1:NB
            temp = I_ns_nb{k,b};
            I_ns_n = I_ns_n + temp;
        end
        FIM0 = 2/c.sigma^2*real(I_ns_n);
        FIM = FIM + FIM0;
    end
    FIM = J*FIM*J.';

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
    if c.syn_type == "Syn"
    	CRLB = FIM(1:end-1,1:end-1)^(-1);
    else
        CRLB = FIM^(-1);
    end

%     sqrt(diag(CRLB(1:3, 1:3)))

    
    
end