function cf = get_CRLB_Ideal(cf)
    snr = cf.snr;
    P_vec = cf.P_vec;
    N0_f_all = 1./db2pow(snr)./db2pow(P_vec);
    % Derivative: exp(1j*phi)*Ast
    FIM_M = zeros(4,4);

    % exp(-1j*2*pi*(0:cf.K-1)'*cf.delta_f*cf.d0/cf.c); % {K x 1}
    H = cf.gain*cf.AstRx*cf.AstTx.';
    D_Adelay_d = cf.Adelay.*(-1j*2*pi*(0:cf.K-1)'*cf.delta_f/cf.c);
    D_theta_H = cf.gain*exp(1j*pi*cf.D_Rx*sind(cf.AOA0))*1j*pi.*cf.D_Rx*cosd(cf.AOA0)*cf.AstTx.';
    D_rho_H = cf.gain/cf.rho*exp(1j*pi*cf.D_Rx*sind(cf.AOA0))*cf.AstTx.';
    D_xi_H = 1j*cf.gain*cf.AstRx*cf.AstTx.';

    for g = 1:cf.G
        Du_theta = cf.WRx(:, g).'*cf.C_Rx*D_theta_H*cf.C_Tx*cf.WTx(:, g)*(cf.Adelay.*cf.s);
        Du_d = cf.WRx(:, g).'*cf.C_Rx*H*cf.C_Tx*cf.WTx(:, g)*(D_Adelay_d.*cf.s);
        Du_rho = cf.WRx(:, g).'*cf.C_Rx*D_rho_H*cf.C_Tx*cf.WTx(:, g)*(cf.Adelay.*cf.s); 
        Du_xi = cf.WRx(:, g).'*cf.C_Rx*D_xi_H*cf.C_Tx*cf.WTx(:, g)*(cf.Adelay.*cf.s);
        Du = [Du_theta, Du_d, Du_rho, Du_xi];
        FIM_M = FIM_M + real(Du'*Du);
    end

    p.tBU_loc = get_dir_from_angle([cf.AOA0 0]);
    p.dBU = cf.d0;

    p.RotB = eye(3);
    D_PU_theta = (p.tBU_loc(1)*p.RotB(:,2) - p.tBU_loc(2)*p.RotB(:,1))/(p.tBU_loc(1)^2 + p.tBU_loc(2)^2)/p.dBU;
    D_PU_d = [cf.PU0; 0]/p.dBU;
    J = eye(4);
    J(1:2, 1:2) = [D_PU_theta(1:2) D_PU_d(1:2)];
    % real(Du_theta*Du_theta')
    % FIM0 = real(Du'*Du);
    % EFIM = FIM0(1,1)-FIM0(1,2:end)*FIM0(2:end,2:end)^-1*FIM0(2:end,1);
    CRLB_M = (FIM_M)^-1;

    cf.AEB_Ideal = rad2deg(sqrt(CRLB_M(1,1)*N0_f_all/2));
    cf.DEB_Ideal = sqrt(CRLB_M(2,2)*N0_f_all/2);

    % FIM = 2./N0_f_all*EFIM;
    FIM_S = J*FIM_M*J.';
    CRB_S = (FIM_S^-1);
    cf.PEB_Ideal = sqrt(trace(CRB_S(1:2, 1:2))*N0_f_all/2);

end