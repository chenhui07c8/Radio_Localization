function CRLB = get_crlb(p)
    % Jacobian matrix
    p.J = get_jacobian_matrix(p);

    % Rotation over element position
    p.D_PhiM_pm = get_rotation_over_position(p);

    J = p.J;
    D_PhiM_pm = p.D_PhiM_pm;
    
    N = p.N;
    NB = p.NB;
    NM = p.NM;
    NR = p.NR;
    B = p.B;
    M = p.M;
    R = p.R;
    P = p.P;
    c = p.c;
    Xmn = p.Xmn;
    fn = p.fn;
    thetaBM = p.thetaBM;
    phiBM = p.phiBM;
    dBM = p.dBM;
    PM = p.PM;
    PB = p.PB;
    tBM = p.tBM;
    rhoBM = p.rhoBM;
    Omega = p.Omega;
    thetaRM = p.thetaRM;
    phiRM = p.phiRM;
    dRM = p.dRM;
    PR = p.PR;
    tRM = p.tRM;
    rhoBRM = p.rhoBRM;
    BW = p.BW;
    KR = p.KR;
    I_ns_nb = cell(N,NB);

    AeqBM = p.AeqBM;
    D_phi_AeqBM = p.D_phi_AeqBM;
    D_theta_AeqBM = p.D_theta_AeqBM;
    
    AeqBR = p.AeqBR;
    
    AeqMB = p.AeqMB;
    D_phi_AeqMB = p.D_phi_AeqMB;
    D_theta_AeqMB = p.D_theta_AeqMB;
    D_alpha_AeqMB = p.D_alpha_AeqMB;
    D_beta_AeqMB = p.D_beta_AeqMB;
    D_gamma_AeqMB = p.D_gamma_AeqMB;
    

    AeqMR = p.AeqMR;
    D_phi_AeqMR = p.D_phi_AeqMR;
    D_theta_AeqMR = p.D_theta_AeqMR;
    D_alpha_AeqMR = p.D_alpha_AeqMR;
    D_beta_AeqMR = p.D_beta_AeqMR;
    D_gamma_AeqMR = p.D_gamma_AeqMR;
    
    for n = 1:N
        % loop BS antennas
        for b = 1:NB
            % taubm = vecnorm(B(:, b)-M,2,1);
            ubm = zeros(NM, 1);
            D_thetaBM_taubm = zeros(NM, 1);
            D_phiBM_taubm = zeros(NM, 1);
            D_tauBM_taubm = zeros(NM, 1);
            ubrm = zeros(NM, NR);
            D_thetaRM_taurm = zeros(NM, NR);
            D_phiRM_taurm = zeros(NM, NR);
            D_tauRM_taurm = zeros(NM, NR);

            % Phi
            part1 = 0;
            part2 = 0;
            for m = 1:NM
                % BM *******
                taubm = norm(B(:,b) - M(:,m))/c;
                dbm = norm(B(:,b)-M(:,m));
                ubm(m) = Xmn(n, m)*exp(-1j*2*pi*fn(n)*taubm);
    %             dbm = taubmMat(b, m)*c;
    %             ubm(m) = Xmn(n, m)*exp(-1j*2*pi*fn(n)*(taubmMat(b, m)));
                % theta
                t_temp = [sin(thetaBM)*cos(phiBM); sin(thetaBM)*sin(phiBM); -cos(thetaBM)];
                D_thetaBM_taubm(m) = dBM/c/dbm*((M(:,m)-PM-B(:,b)+PB)'*t_temp);
                % phi
                t_temp = [cos(thetaBM)*sin(phiBM); -cos(thetaBM)*cos(phiBM); 0];
                D_phiBM_taubm(m) = dBM/c/dbm*((M(:,m)-PM-B(:,b)+PB)'*t_temp);
                % tau
                Gbm2 = (M(:,m)-PM - B(:,b)+PB)'*tBM;
                D_tauBM_taubm(m) = 1/dbm*(dBM - Gbm2);

                % PhiM part 1
                D_PhiM_taubm = 1/c/dbm*((M(:,m)-PM-B(:,b)+PB)'*D_PhiM_pm{m} - dBM*(tBM'*D_PhiM_pm{m}));
                part1 = part1 + rhoBM*ubm(m)*D_PhiM_taubm;

                % BRM ********
                for r = 1:NR
                    taubr = norm(B(:,b) - R(:,r))/c;
                    taurm = norm(R(:,r) - M(:,m))/c;
                    drm = norm(R(:,r)-M(:,m));
                    ubrm(m,r) = Xmn(n, m)*Omega(r)*exp(-1j*2*pi*fn(n)*(taubr + taurm));
    %                 drm = taurmMat(r, m)*c;
    %                 ubrm(m,r) = Xmn(n, m)*Omega(n, r)*exp(-1j*2*pi*fn(n)*(taubrMat(b, r) + taurmMat(r, m)));
                    % theta
                    t_temp = [sin(thetaRM)*cos(phiRM); sin(thetaRM)*sin(phiRM); -cos(thetaRM)];
                    D_thetaRM_taurm(m,r) = dRM/c/drm*((M(:,m)-PM-R(:,r)+PR)'*t_temp);
                    % phi
                    t_temp = [cos(thetaRM)*sin(phiRM); -cos(thetaRM)*cos(phiRM); 0];
                    D_phiRM_taurm(m,r) = dRM/c/drm*((M(:,m)-PM-R(:,r)+PR)'*t_temp);
                    % tau
                    Grm2 = (M(:,m)-PM - R(:,r)+PR)'*tRM;
                    D_tauRM_taurm(m,r) = 1/drm*(dRM - Grm2);

                    % PhiM part2
                    D_PhiM_taurm = 1/c/drm*((M(:,m)-PM-R(:,r)+PR)'*D_PhiM_pm{m} - dRM*(tRM'*D_PhiM_pm{m}));
                    part2 = part2 + rhoBRM*ubrm(m,r)*D_PhiM_taurm;
                end
            end
            ubBM = sqrt(P)*rhoBM*sum(ubm);
            ubBRM = sqrt(P)*rhoBRM*sum(ubrm(:));

            % two-stage localization
            % BM
            D_rhoBM_ub = ubBM/rhoBM;
            D_thetaBM_ub = -1j*2*pi*fn(n)*sqrt(P)*rhoBM*sum(ubm.*D_thetaBM_taubm);
            D_phiBM_ub = - 1j*2*pi*fn(n)*sqrt(P)*rhoBM*sum(ubm.*D_phiBM_taubm);
            D_tauBM_ub = - 1j*2*pi*fn(n)*sqrt(P)*rhoBM*sum(ubm.*D_tauBM_taubm);

            % BRM
            D_rhoBRM_ub = ubBRM/rhoBRM;
            D_thetaRM_ub = -1j*2*pi*fn(n)*sqrt(P)*rhoBRM*sum(ubrm.*D_thetaRM_taurm, 'all');
            D_phiRM_ub   = -1j*2*pi*fn(n)*sqrt(P)*rhoBRM*sum(ubrm.*D_phiRM_taurm, 'all');
            D_tauRM_ub   = -1j*2*pi*fn(n)*sqrt(P)*rhoBRM*sum(ubrm.*D_tauRM_taurm, 'all');

            D_PhiM_ub_1 = -1j*2*pi*fn(n)*sqrt(P)*part1;
            D_PhiM_ub_2 = -1j*2*pi*fn(n)*sqrt(P)*part2;
            D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;

            
            % extra calculation with AOSA
            
            % PWM for SA AOA/AOD
            % BM channel
            D_phiBM_ub = D_phiBM_ub*AeqBM(b, n)*AeqMB(b, n) + ubBM*D_phi_AeqBM(b, n)*AeqMB(b, n) + ubBM*AeqBM(b, n)*D_phi_AeqMB(b, n);
            D_thetaBM_ub = D_thetaBM_ub*AeqBM(b, n)*AeqMB(b, n) + ubBM*D_theta_AeqBM(b, n)*AeqMB(b, n) + ubBM*AeqBM(b, n)*D_theta_AeqMB(b, n);
            D_PhiM_AeqMB = [D_alpha_AeqMB(b, n), D_beta_AeqMB(b, n), D_gamma_AeqMB(b, n)];
            D_PhiM_ub_1 = AeqBM(b, n)*(D_PhiM_ub_1*AeqMB(b, n) + ubBM*D_PhiM_AeqMB);
            % BRM channel
            D_phiRM_ub = AeqBR(b, n)*(D_phiRM_ub*AeqMR(b,n) + ubBRM*D_phi_AeqMR(b, n));
            D_thetaRM_ub = AeqBR(b, n)*(D_thetaRM_ub*AeqMR(b,n) + ubBRM*D_theta_AeqMR(b, n));
            D_PhiM_AeqMR = [D_alpha_AeqMR(b, n), D_beta_AeqMR(b, n), D_gamma_AeqMR(b, n)];
            D_PhiM_ub_2 = AeqBR(b, n)*(D_PhiM_ub_2*AeqMR(b, n) + ubBRM*D_PhiM_AeqMR);

            D_PhiM_ub = D_PhiM_ub_1 + D_PhiM_ub_2;
            
            % SWM for SA AOA/AOD
            
            
            
            
            % done with calculations 
            D_ub = zeros(1,11);
            D_ub(1:4) = [D_rhoBM_ub, D_thetaBM_ub, D_phiBM_ub, D_tauBM_ub];
            D_ub(5:8) = [D_rhoBRM_ub, D_thetaRM_ub, D_phiRM_ub, D_tauRM_ub];
            D_ub(9:11) = D_PhiM_ub;

            I_ns_nb{n,b} = (D_ub'*D_ub);
        end
    end


    FIM = zeros(6,6);
    for n = 1:N
        I_ns_n = zeros(11,11);
        for b = 1:NB
            I_ns_n = I_ns_n + I_ns_nb{n,b};
        end
        FIM0 = 2/p.sigma^2*real(I_ns_n);
        FIM = FIM + J*FIM0*J.';
    end

    % get crlb
    CRLB = FIM^(-1);


end