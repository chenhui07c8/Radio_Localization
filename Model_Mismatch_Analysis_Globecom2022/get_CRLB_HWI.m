function cp = get_CRLB_HWI(cp)
    P_vec = cp.P_vec;
    snr = cp.snr;
    F = fft(eye(cp.K))/sqrt(cp.K);
    FH = ifft(eye(cp.K))*sqrt(cp.K);
%     N0_p_all = zeros(1, length(P_vec));
%     X_cell = cell(1, length(P_vec));
    CRB_p_angle = zeros(1, length(P_vec));
    CRB_p_delay = zeros(1, length(P_vec));
    CRB_p_pos = zeros(1, length(P_vec));

    for P_i = 1:length(P_vec)
        cp.P = db2pow(cp.P_vec(P_i));
        cp.PU = cp.PU0;
        cp.AOA = cp.AOA0;
        cp.AOD = cp.AOD0;
        cp.d = cp.d0;

        cp = get_Rx_symbol_HWI(cp);
        %     X_cell{P_i} = cp.X;
        N0_p_all = 1./db2pow(snr);
        %     N0_p_all(P_i) = 1./db2pow(snr)./norm(cp.X(:, 1)).^2;

        % get derivative w.r.t the channel
        FIM_M = zeros(4,4);
        H = cp.gain*cp.AstRx*cp.AstTx.';
        D_Adelay_D = cp.Adelay.*(-1j*2*pi*(0:cp.K-1)'*cp.delta_f/cp.c);
        D_theta_H = cp.gain*exp(1j*pi*cp.D_Rx*sind(cp.AOA0))*1j*pi.*cp.D_Rx*cosd(cp.AOA0)*cp.AstTx.';
        D_rho_H = cp.gain/cp.rho*exp(1j*pi*cp.D_Rx*sind(cp.AOA0))*cp.AstTx.';
        D_xi_H = 1j*cp.gain*cp.AstRx*cp.AstTx.';

        % get derivative w.r.t the received symbol
        for g = 1:cp.G
            D_tildeu_theta = (cp.vec_CFO_Rx.*cp.vec_PN_Rx(:, g)).*...
                ( FH* (cp.WRx(:, g).'*cp.C_Rx*D_theta_H*cp.C_Tx*cp.WTx(:, g)*(cp.Adelay.*cp.X)) );
            D_tildeu_d = (cp.vec_CFO_Rx.*cp.vec_PN_Rx(:, g)).*...
                ( FH* (cp.WRx(:, g).'*cp.C_Rx*H*cp.C_Tx*cp.WTx(:, g)*(D_Adelay_D.*cp.X)) );
            D_tildeu_rho = (cp.vec_CFO_Rx.*cp.vec_PN_Rx(:, g)).*...
                ( FH* (cp.WRx(:, g).'*cp.C_Rx*D_rho_H*cp.C_Tx*cp.WTx(:, g)*(cp.Adelay.*cp.X)) ); 
            D_tildeu_xi = (cp.vec_CFO_Rx.*cp.vec_PN_Rx(:, g)).*...
                ( FH* (cp.WRx(:, g).'*cp.C_Rx*D_xi_H*cp.C_Tx*cp.WTx(:, g)*(cp.Adelay.*cp.X)) );


            Du_theta=   F*cp.IQI_alpha_Rx*D_tildeu_theta + F*cp.IQI_beta_Rx*conj(D_tildeu_theta);
            Du_d    =   F*cp.IQI_alpha_Rx*D_tildeu_d + F*cp.IQI_beta_Rx*conj(D_tildeu_d);
            Du_rho  =   F*cp.IQI_alpha_Rx*D_tildeu_rho + F*cp.IQI_beta_Rx*conj(D_tildeu_rho);
            Du_xi   =   F*cp.IQI_alpha_Rx*D_tildeu_xi + F*cp.IQI_beta_Rx*conj(D_tildeu_xi);
            D_u = [Du_theta, Du_d, Du_rho, Du_xi];
            FIM_M = FIM_M + real(D_u'*D_u);
        end

        % get Jacobian matrix from Measurement vector to State vector
        p.tBU_loc = get_dir_from_angle([cp.AOA0 0]);
        p.dBU = cp.d0;
        p.RotB = eye(3);
        D_PU_theta = (p.tBU_loc(1)*p.RotB(:,2) - p.tBU_loc(2)*p.RotB(:,1))/(p.tBU_loc(1)^2 + p.tBU_loc(2)^2)/p.dBU;
        D_PU_d = [cp.PU0; 0]/p.dBU;
        J = eye(4);
        J(1:2, 1:2) = [D_PU_theta(1:2) D_PU_d(1:2)];
        
        % CRLB of the Measurement Vector and the State Vector
        % real(Du_theta*Du_theta')
        % FIM0 = real(Du'*Du);
        % EFIM = FIM0(1,1)-FIM0(1,2:end)*FIM0(2:end,2:end)^-1*FIM0(2:end,1);
        CRLB_M = (FIM_M)^-1;
        CRB_p_angle(P_i) = rad2deg(sqrt(CRLB_M(1,1)*N0_p_all/2));
        CRB_p_delay(P_i) = sqrt(CRLB_M(2,2)*N0_p_all/2);

        % FIM = 2./N0_f_all*EFIM;
        FIM_S = J*FIM_M*J.';
        CRB_S = (FIM_S^-1);
        CRB_p_pos(P_i) = sqrt(trace(CRB_S(1:2, 1:2))*N0_p_all/2);
        
        cp.AEB_HWI = CRB_p_angle;
        cp.DEB_HWI = CRB_p_delay;
        cp.PEB_HWI = CRB_p_pos;
        
    end

    
end