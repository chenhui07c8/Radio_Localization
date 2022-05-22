

function FIM0 = get_FIM_PWM(c)

%     c.D_PhiM_pm = get_rotation_over_position(c);
%     D_PhiM_pm = c.D_PhiM_pm;
    if c.AoSA_type == "SPP"
        c = get_D_sv_Aeq_SPP(c);
    elseif c.AoSA_type == "SSP"
        c = get_D_sv_Aeq_SSP(c);
    end
    
    B = c.B;
    M = c.M;
    K = c.K;
    NB = c.NB;
    NM = c.NM;
    NR = c.NR;

    thetaBM = c.thetaBM;
    phiBM = c.phiBM;
    dBM = c.dBM;
    PM = c.PM;
    PB = c.PB;
    tBM = c.tBM;
%     rhoBM = c.rhoBM;

    if(c.LR > 0)
        R = c.R;
        thetaRM = c.thetaRM;
        phiRM = c.phiRM;
        dRM = c.dRM;
        PR = c.PR;
        tRM = c.tRM;
    end
%     rhoBRM = c.rhoBRM;

    P = c.P;
    fk = c.fk;
    I_ns_k = cell(K);

    % derivative matrix of BM channel
    [D_phiBM_tBM, D_thetaBM_tBM] = get_D_Phi_t(c.phiBM, c.thetaBM);
    % BM
    OmegaPsiBM = c.B0'*c.RotB'*c.tBM;
    AstBM = exp(1j*2*pi./c.lambdak.*(OmegaPsiBM));
    D_phiBM_AstBM = AstBM.*(1j*2*pi./c.lambdak.*(c.B0'*c.RotB'*D_phiBM_tBM));
    D_thetaBM_AstBM = AstBM.*(1j*2*pi./c.lambdak.*(c.B0'*c.RotB'*D_thetaBM_tBM));
    
    rhoBM = c.rhoBM;

    
    AeqBM = c.AeqBM;
    AeqMB = c.AeqMB;

    D_phiBM_AeqBM = c.D_phiBM_AeqBM;
    D_phiBM_AeqMB = c.D_phiBM_AeqMB;
    D_thetaBM_AeqBM = c.D_thetaBM_AeqBM;
    D_thetaBM_AeqMB = c.D_thetaBM_AeqMB;
    % calculate FIM
    for k = 1:K
        D_beta_ub = zeros(c.NB, 3);
        % BM channel 
        if (c.LB > 0)
            ubm = c.XiL(:,:,k);
            ubBM = c.uBM(:, k);
            D_rhoBM_ub = ubBM/rhoBM;
            D_xiBM_ub = 1j*2*pi/c.lambdac*ubBM;
            D_beta_ub(:, 1) = -1j*2*pi/c.lambdak(k)*(ubBM);

            % AOSA structure
            % BM channel
            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*D_phiBM_AstBM(:,k).*transpose(AeqMB(:,:,k)) + ...
                D_phiBM_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                AeqBM(:,:,k).*ubm.*transpose(D_phiBM_AeqMB(:,:,k));
            D_phiBM_ub = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);

            temp = AeqBM(:,:,k)*-1j*2*pi*fk(k).*ubm.*D_thetaBM_AstBM(:,k).*transpose(AeqMB(:,:,k)) + ...
                    D_thetaBM_AeqBM(:,:,k).*ubm.*transpose(AeqMB(:,:,k)) + ...
                    AeqBM(:,:,k).*ubm.*transpose(D_thetaBM_AeqMB(:,:,k));
            D_thetaBM_ub = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);

            temp = AeqBM(:,:,k).*(-1j*2*pi*c.fdk(k)*ubm).*transpose(AeqMB(:,:,k));
            D_tauBM_ub = sqrt(P)*c.alphaL(k)*temp*c.Xmk(:,k);
        end
        
        D_beta_ub = sum(D_beta_ub, 2);
        
        D_ub = zeros(NB, c.LB*5 + 3 +1);
        if (c.LB > 0)
            for l = 1:c.LB
            	D_ub(:, (1:5) + 5*(l-1)) = [D_rhoBM_ub(:, l),  D_xiBM_ub(:, l),  D_phiBM_ub(:, l),...
                    D_thetaBM_ub(:, l), D_tauBM_ub(:, l)];
            end
        end
        D_ub(:, end) = D_beta_ub;

%         if (c.LR > 0)
%             for l = 1:c.LR
%             	D_ub(:, 5*c.LB + (1:5) + 5*(l-1)) = [D_rhoBRM_ub(:, l),  D_xiBRM_ub(:, l),  D_phiRM_ub(:, l),...
%                     D_thetaRM_ub(:, l), D_tauRM_ub(:, l)];
%             end
%         end
%         D_ub(:, end-3:end-1) = D_PhiM_ub;
%         if (c.LC > 0)
%             % each cluster has 6 unknowns in 2D (rho, xi, phiBC, tauBC, phiCM, tauCM)
%             for l = 1:c.LC
%                 D_ub(:, (5*c.LB + 5*c.LR +(1:8) + 8*(l-1))) = [D_rhoBCM_ub(:,l),  D_xiBCM_ub(:,l),  D_phiBC_ub(:,l), D_thetaBC_ub(:,l), ...
%                     D_tauBC_ub(:,l), D_phiCM_ub(:,l), D_thetaCM_ub(:,l), D_tauCM_ub(:,l)];
%             end
%         end
        

        I_ns_k{k} = (D_ub'*D_ub);  
    end

    
    FIM0 = zeros(size(D_ub,2), size(D_ub,2));
    for k = 1:K
        I_ns_n = zeros(size(FIM0));
        temp = I_ns_k{k};
        I_ns_n = I_ns_n + temp;
        FIM0 = FIM0 + 2/c.sigma^2*real(I_ns_n);
    end





end
