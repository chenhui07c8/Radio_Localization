function [AEB_LB, DEB_LB, PEB_LB, eta_opt] = get_CRLB_PWM_LB(cp, cf)

    P_vec = cp.P;
    snr = cp.snr;
    MCRB_angle = zeros(size(P_vec));
    MCRB_delay = zeros(size(P_vec));

    AEB_LB = zeros(size(P_vec));
    DEB_LB = zeros(size(P_vec));
    PEB_LB = zeros(size(P_vec));
    
    PU0 = cp.PU;
    eta = [deg2rad(cp.AOA) cp.d cp.rho mod(cp.xi, 2*pi)].';

    for P_i = 1:length(P_vec)
        % P_i = length(P_vec);
%         disp(num2str(P_i) + " of " + num2str(length(P_vec)));
        P_cur = (P_vec(P_i));
        
        cp.P = P_cur;
        cp.PU = cp.PU;
        cp.AOA = atan2d(cp.PU(2), cp.PU(1));    % angle of departure
        cp.d = norm(cp.PU);
        cp = get_Rx_symbol(cp);
        yp = cp.u;
        
        cf.P = P_cur;
        cf = get_Rx_symbol(cf);
        % yf = cf.u; norm(yf-yp, 'fro')
        [PU_opt, gain_opt] = alg_mle_toy_POS(yp, cf);
        rho_opt = abs(gain_opt);
        xi_opt = phase(gain_opt);
        theta_opt = atan2d(PU_opt(2), PU_opt(1));
        d_opt = norm(PU_opt);

        AstTx = exp(1j*pi*cf.D_Tx*sind(0));  % steering vector at the Rx
        AstRx = exp(1j*pi*cf.D_Rx*sind(theta_opt));  % steering vector at the Rx
        H_opt = gain_opt*AstRx*AstTx.';
%         Adelay = exp(-1j*2*pi*(0:cf.K-1)*cf.delta_f*d_opt/cf.c); % {K x 1}
        Adelay = exp(-1j*2*pi./cf.lambdak*d_opt);
        yf_opt = zeros(cf.M_Rx, cf.K, cf.G);
        for g = 1:cf.G
            yf_opt(:, :, g) = cf.WRx(:, g).'*H_opt*cf.WTx(:, g)*(Adelay.*cf.X); % cf.G*cf.K
        end
        % norm(yf_opt - yp, 'fro')
        
        % MCLB
        % variables [theta, d, rho, xi]
        D_mu_cell = cell(1, 4); % theta, d, rho, xi
        % each cell has the size of K x G
        D_theta_H = gain_opt*exp(1j*pi*cf.D_Rx*sind(theta_opt))*1j*pi.*cf.D_Rx*cosd(theta_opt)*AstTx.';
        D_Adelay_d = Adelay.*(-1j*2*pi./cf.lambdak);
        D_rho_H = exp(1j*xi_opt)*exp(1j*pi*cf.D_Rx*sind(theta_opt))*AstTx.';
        D_xi_H = 1j*gain_opt*exp(1j*pi*cf.D_Rx*sind(theta_opt))*AstTx.';
        for g = 1:cf.G
            D_mu_cell{1}(:, g) = reshape(cf.WRx(:, g).'*D_theta_H*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            D_mu_cell{2}(:, g) = reshape(cf.WRx(:, g).'*H_opt*cf.WTx(:, g)*(D_Adelay_d.*cf.X), [cf.K*cf.M_Rx 1]);
            D_mu_cell{3}(:, g) = reshape(cf.WRx(:, g).'*D_rho_H*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            D_mu_cell{4}(:, g) = reshape(cf.WRx(:, g).'*D_xi_H*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);

            D2_mu_cell = cell(4, 4);
            % theta: theta, d, rho, xi
            D2_mu_cell{1,1}(:, g) = reshape(cf.WRx(:, g).'*(D_theta_H.*(1j*pi*cf.D_Rx*cosd(theta_opt))...
                    - gain_opt*exp(1j*pi*cf.D_Rx*sind(theta_opt)).*(1j*pi*cf.D_Rx*sind(theta_opt))*AstTx.')...
                    *cf.WTx(:, g)*(cf.Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{1,2}(:, g) = reshape(cf.WRx(:, g).'*D_theta_H*cf.WTx(:, g)*(D_Adelay_d.*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{1,3}(:, g) = reshape(cf.WRx(:, g).'*(exp(1j*xi_opt)*exp(1j*pi*cf.D_Rx*sind(theta_opt))...
                    .*(1j*pi*cf.D_Rx*cosd(theta_opt)))*AstTx.'*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{1,4}(:, g) = reshape(cf.WRx(:, g).'*(1j*gain_opt*exp(1j*pi*cf.D_Rx*sind(theta_opt))...
                    .*(1j*pi*cf.D_Rx*cosd(theta_opt)))*AstTx.'*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            % d: theta, d, rho, xi
            D2_mu_cell{2,2}(:, g) = reshape(cf.WRx(:, g).'*H_opt*cf.WTx(:, g)...
                *(D_Adelay_d.*(-1j*2*pi./cf.lambdak).*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{2,3}(:, g) = reshape(cf.WRx(:, g).'*D_rho_H*cf.WTx(:, g)*(D_Adelay_d.*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{2,4}(:, g) = reshape(cf.WRx(:, g).'*D_xi_H*cf.WTx(:, g)*(D_Adelay_d.*cf.X), [cf.K*cf.M_Rx 1]);
            % rho: rho, xi
            D2_mu_cell{3,3}(:, g) = zeros(cf.K*cf.M_Rx, 1);
            D2_mu_cell{3,4}(:, g) = reshape(cf.WRx(:, g).'*...
                    (1j*exp(1j*xi_opt)*exp(1j*pi*cf.D_Rx*sind(theta_opt)))*AstTx.'*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{4,4}(:, g) = reshape(cf.WRx(:, g).'*...
                    (-1*gain_opt*exp(1j*pi*cf.D_Rx*sind(theta_opt)))*AstTx.'*cf.WTx(:, g)*(Adelay.*cf.X), [cf.K*cf.M_Rx 1]);
            D2_mu_cell{2,1}(:, g) = D2_mu_cell{1,2}(:, g);
            D2_mu_cell{3,1}(:, g) = D2_mu_cell{1,3}(:, g);
            D2_mu_cell{3,2}(:, g) = D2_mu_cell{2,3}(:, g);
            D2_mu_cell{4,1}(:, g) = D2_mu_cell{1,4}(:, g);
            D2_mu_cell{4,2}(:, g) = D2_mu_cell{2,4}(:, g);
            D2_mu_cell{4,3}(:, g) = D2_mu_cell{3,4}(:, g);
        end
        epsilon = yp - yf_opt;
%         norm(yp-yf_opt, 'fro')
% yp()
        eta_opt = [deg2rad(theta_opt) d_opt rho_opt mod(xi_opt, 2*pi)].';

        N0 = 1/db2pow(snr);
%     N0 = 1
%         N0 = 1/db2pow(snr)/db2pow(P_vec(P_i));
        A = zeros(4,4);
        B = zeros(4,4);
        for i = 1:4
            for j = 1:4
                A(i,j) = 2/N0*real(epsilon(:)'*D2_mu_cell{i,j}(:) - D_mu_cell{i}(:)'*D_mu_cell{j}(:));
                B(i,j) = 2/N0*(2/N0*real(epsilon(:)'*D_mu_cell{i}(:))*real(epsilon(:)'*D_mu_cell{j}(:)) + ...
                        real(D_mu_cell{i}(:)'*D_mu_cell{j}(:)));
            end
        end
        
        % standard version
%         M_CRB = A^-1*B*A^-1;
        
        % scaled amplitude for better calculation
        amp_scale = 1000;
        A2 = A;
        A2(3, :) = A2(3, :)/amp_scale;
        A2(:, 3) = A2(:, 3)/amp_scale;
        B2 = B;
        B2(3, :) = B2(3, :)/amp_scale;
        B2(:, 3) = B2(:, 3)/amp_scale;
        M_CRB = A2^-1*B2*A2^-1;

    %     M_CRB = -A^-1;    % if no mismatch
        MCRB_angle(P_i) = rad2deg(sqrt(trace(M_CRB(1,1))));
        MCRB_delay(P_i) = (sqrt(trace(M_CRB(2,2))));

        L_CRB = M_CRB + (eta - eta_opt)*(eta - eta_opt)';
        AEB_LB(P_i) = rad2deg(sqrt(trace(L_CRB(1,1))));
        DEB_LB(P_i) = (sqrt(trace(L_CRB(2,2))));

        cf.tBU_loc = get_dir_from_angle([rad2deg(eta(1)) 0]');
        cf.dBU = eta(2);
        cf.RotB = eye(3);
        D_PU_theta = (cf.tBU_loc(1)*cf.RotB(:,2) - cf.tBU_loc(2)*cf.RotB(:,1))/(cf.tBU_loc(1)^2 + cf.tBU_loc(2)^2)/cf.dBU;
        D_PU_d = [PU0; 0]/cf.dBU;
        J = eye(4);
        J(1:2, 1:2) = [D_PU_theta(1:2) D_PU_d(1:2)];
%         CRLB_S = (J*L_CRB^-1*J.')^-1;
        CRLB_S = (J.')^-1*L_CRB*J^-1;
        PEB_LB(P_i) = (sqrt(trace(CRLB_S(1:2,1:2))));

        
%         J = [D_PU_theta(1:2) D_PU_d(1:2)];
%         E_LB = L_CRB(1:2,1:2);
%         CRLB_S = (J*E_LB^-1*J.')^-1;
%         PEB_LB(P_i) = (sqrt(trace(CRLB_S(1:2,1:2))));
        
        
    %     eta = [deg2rad(cp.AOA0) cp.d0 cp.rho cp.xi].';
    %     eta_opt = [deg2rad(theta_opt) d_opt rho_opt xi_opt].';
    %     LS_CRB = J*M_CRB^-1*J.';
    %     CRLB_S = (LS_CRB^-1);
    %     LB_pos(snr_i) = (sqrt(trace(CRLB_S(1:2,1:2))));

    end
end