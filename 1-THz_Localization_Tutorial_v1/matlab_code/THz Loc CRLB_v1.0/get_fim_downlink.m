function [FIM0] = get_fim_downlink(c)
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
    
    % get derivatives of matrix
    c = get_D_measure_from_matrix(c);
    D_thetaBM_taubm = c.D_thetaBM_taubm;
    D_phiBM_taubm = c.D_phiBM_taubm;
    D_tauBM_taubm = c.D_tauBM_taubm;
    D_PhiM_ubBM = c.D_PhiM_ubBM;
    
    D_thetaRM_taurm = c.D_thetaRM_taurm;
    D_phiRM_taurm = c.D_phiRM_taurm;
    D_tauRM_taurm = c.D_tauRM_taurm;
    D_PhiM_ubBRM = c.D_PhiM_ubBRM;

    
    K = c.K;
    NB = c.NB;
    NM = c.NM;
    NR = c.NR;
    P = c.P;
    fk = c.fk;
    rhoBM = c.rhoBM;
    rhoBRM = c.rhoBRM;
    I_ns_k = cell(K);

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
    
    
    % calculate FIM
    for k = 1:K
        ubm = c.XiL(:,:,k);
        ubr = c.XiBR(:,:,k);
        urm = c.XiRM(:,:,k);
        
        ubBM = c.uMB(:, k);
        ubBRM = c.uMRB(:, k);
        D_rhoBM_ub = ubBM/rhoBM;
        D_xiBM_ub = 1j*2*pi/c.lambdac*ubBM;
        D_rhoBRM_ub = ubBRM/rhoBRM;
        D_xiBRM_ub = 1j*2*pi/c.lambdac*ubBRM;
        D_beta_ub = -1j*2*pi/c.lambdak(k)*(ubBM + ubBRM);
        % AOSA structure
        
        % BM channel (NM x 1)
        temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*D_phiBM_taubm(:,:,k).*transpose(AeqMB(:,:,k)) + ...
            D_phiBM_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
            AeqBM(:,:,k).*ubm.*transpose(D_phiBM_AeqMB(:,:,k));
        D_phiBM_ub = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);
            
        temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*D_thetaBM_taubm(:,:,k).*transpose(AeqMB(:,:,k)) + ...
                D_thetaBM_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubm.*transpose(D_thetaBM_AeqMB(:,:,k));
        D_thetaBM_ub = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);
        
        temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*D_tauBM_taubm(:,:,k).*transpose(AeqMB(:,:,k)) + ...
                D_tauBM_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubm.*transpose(D_tauBM_AeqMB(:,:,k));
        D_tauBM_ub = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);
        % orientation ub1
        temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*reshape(D_PhiM_ubBM(:,:,k,1),[NB NM]).*transpose(AeqMB(:,:,k)) + ...
                D_alpha_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubm.*transpose(D_alpha_AeqMB(:,:,k));
        D_alpha_ub_1 = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

        temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*reshape(D_PhiM_ubBM(:,:,k,2),[NB NM]).*transpose(AeqMB(:,:,k)) + ...
                D_beta_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubm.*transpose(D_beta_AeqMB(:,:,k));
        D_beta_ub_1 = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);

        temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*reshape(D_PhiM_ubBM(:,:,k,3),[NB NM]).*transpose(AeqMB(:,:,k)) + ...
                D_gamma_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubm.*transpose(D_gamma_AeqMB(:,:,k));
        D_gamma_ub_1 = sqrt(P)*c.alphaL(k)*transpose(temp)*c.Xbk(:,k);
        D_PhiM_ub_1 = [D_alpha_ub_1 D_beta_ub_1 D_gamma_ub_1];
        
        
        % BRM channel
        temp = AeqBR(:,:,k).*ubr*diag(c.Omega)* ...
            transpose(-1j*2*pi*fk(k)*urm.*D_phiRM_taurm(:,:,k).*AeqMR(:,:,k) + urm.*D_phiRM_AeqMR(:,:,k));
        D_phiRM_ub = sqrt(P)*c.alphaR(k)*transpose(temp)*c.Xbk(:,k);
        
        temp = AeqBR(:,:,k).*ubr*diag(c.Omega)*...
            transpose(-1j*2*pi*fk(k)*urm.*D_thetaRM_taurm(:,:,k).*AeqMR(:, :, k) + urm.*D_thetaRM_AeqMR(:,:,k));
        D_thetaRM_ub = sqrt(P)*c.alphaR(k)*transpose(temp)*c.Xbk(:,k);
        
        temp = AeqBR(:,:,k).*ubr*diag(c.Omega)*...
            transpose(-1j*2*pi*fk(k)*urm.*D_tauRM_taurm(:,:,k).*AeqMR(:,:,k) + urm.*D_tauRM_AeqMR(:,:,k));
        D_tauRM_ub = sqrt(P)*c.alphaR(k)*transpose(temp)*c.Xbk(:,k);
        
        % orientation ub2
        temp = AeqBR(:,:,k).*ubr*diag(c.Omega)*...
            transpose(-1j*2*pi*fk(k)*urm.*reshape(D_PhiM_ubBRM(:,:,k,1),[NM NR]).*AeqMR(:,:,k) + urm.*D_alpha_AeqMR(:,:,k));
        D_alpha_ub_2 = sqrt(P)*c.alphaR(k)*transpose(temp)*c.Xbk(:,k);
        temp = AeqBR(:,:,k).*ubr*diag(c.Omega)*...
            transpose(-1j*2*pi*fk(k)*urm.*reshape(D_PhiM_ubBRM(:,:,k,2),[NM NR]).*AeqMR(:,:,k) + urm.*D_beta_AeqMR(:,:,k));
        D_beta_ub_2 = sqrt(P)*c.alphaR(k)*transpose(temp)*c.Xbk(:,k);
        temp = AeqBR(:,:,k).*ubr*diag(c.Omega)*...
            transpose(-1j*2*pi*fk(k)*urm.*reshape(D_PhiM_ubBRM(:,:,k,3),[NM NR]).*AeqMR(:,:,k) + urm.*D_gamma_AeqMR(:,:,k));
        D_gamma_ub_2 = sqrt(P)*c.alphaR(k)*transpose(temp)*c.Xbk(:,k);
        D_PhiM_ub_2 = [D_alpha_ub_2 D_beta_ub_2 D_gamma_ub_2];

        D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;

        if  c.pos_type == "2D"
        % 2D done with calculations 
            D_ub = zeros(NM,10);
            D_ub(:, 1:4) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_tauBM_ub];
            D_ub(:, 5:8) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_tauRM_ub];
            D_ub(:, 9) = D_PhiM_ub(:, 1);
            D_ub(:, 10) = D_beta_ub;
        elseif  c.pos_type == "3D"
        % 3D done with calculations 
            D_ub = zeros(NM,14);
            D_ub(:, 1:5) = [D_rhoBM_ub,  D_xiBM_ub,  D_phiBM_ub, D_thetaBM_ub, D_tauBM_ub];
            D_ub(:, 6:10) = [D_rhoBRM_ub, D_xiBRM_ub, D_phiRM_ub, D_thetaRM_ub, D_tauRM_ub];
            D_ub(:, 11:13) = D_PhiM_ub;
            D_ub(:, 14) = D_beta_ub;
        end
        I_ns_k{k} = (D_ub'*D_ub);  
    end

    
    FIM0 = zeros(size(D_ub,2), size(D_ub,2));
    for k = 1:K
        I_ns_n = zeros(size(FIM0));
        temp = I_ns_k{k};
        I_ns_n = I_ns_n + temp;
        FIM0 = FIM0 + 2/c.sigma^2*real(I_ns_n);
    end
     
%     FIM0 = zeros(length(D_ub), length(D_ub));
%     for k = 1:K
%         I_ns_n = zeros(length(D_ub), length(D_ub));
%         for b = 1:NB
%             temp = I_ns_nb{k,b};
%             I_ns_n = I_ns_n + temp;
%         end
% %         FIM = FIM + J*FIM0*J.';
%         FIM0 = FIM0 + 2/c.sigma^2*real(I_ns_n);
%     end
   
    
    
end