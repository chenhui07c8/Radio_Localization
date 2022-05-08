function cf = get_CRLB_SWM(cf)

    FIM_S = zeros(4,4);

    if cf.NF_SNS == "True"
        % coe_mat = cp.fk./cp.fc.*(norm(cp.PU)./vecnorm(cp.D_Rx_Loc-cp.PU, 2, 1)).';
        dis_temp = vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1).';
        coe_mat = cf.lambdak./cf.lambdac.*(norm(cf.PU)./dis_temp);
        D_coe_x = cf.lambdak./cf.lambdac.*(cf.PU(1)/cf.d*dis_temp - cf.d*(cf.PU(1) - cf.D_Rx_Loc(1,:).')./dis_temp)./(dis_temp.^2);
        D_coe_y = cf.lambdak./cf.lambdac.*(cf.PU(2)/cf.d*dis_temp - cf.d*(cf.PU(2) - cf.D_Rx_Loc(2,:).')./dis_temp)./(dis_temp.^2);
    else
        coe_mat = ones(cf.N, cf.K);
        D_coe_x = zeros(cf.N, cf.K);
        D_coe_y = zeros(cf.N, cf.K);
    end
    
    if cf.NF_SWM == "True"
        dis_vec = vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1).' - norm(cf.PU);
        dvec_mat = exp(-1j*2*pi*dis_vec/cf.lambdac.*cf.bse_coe);
        D_dvec_x = dvec_mat.*(-1j*2*pi/cf.lambdac.*cf.bse_coe).*...
            (((cf.PU(1)-cf.D_Rx_Loc(1,:))./vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1))' - cf.PU(1)/norm(cf.PU));
        D_dvec_y = dvec_mat.*(-1j*2*pi/cf.lambdac.*cf.bse_coe).*...
            (((cf.PU(2)-cf.D_Rx_Loc(2,:))./vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1))' - cf.PU(2)/norm(cf.PU));
    else
        dis_vec = -cf.lambdac/2*(cf.D_Rx-mean(cf.D_Rx))*sind(cf.AOA);
        dvec_mat = exp(-1j*2*pi*dis_vec/cf.lambdac.*cf.bse_coe);
        % sind(cp.AOA) = y/sqrt(x^2 + y^2)
        D_dvec_x = dvec_mat.*(-1j*2*pi/cf.lambdac.*cf.bse_coe).*(-cf.lambdac/2*(cf.D_Rx-mean(cf.D_Rx))*(-cf.PU(2)*cf.PU(1)/norm(cf.PU))/norm(cf.PU)^2);
        D_dvec_y = dvec_mat.*(-1j*2*pi/cf.lambdac.*cf.bse_coe).*(-cf.lambdac/2*(cf.D_Rx-mean(cf.D_Rx))*(norm(cf.PU)-cf.PU(2)^2/norm(cf.PU))/norm(cf.PU)^2);
    end
    
    Adelay = exp(-1j*2*pi*cf.fk*cf.d/cf.c);
    D_Adelay_x = cf.Adelay.*(-1j*2*pi./cf.lambdak*cf.PU(1)/norm(cf.PU));
    D_Adelay_y = cf.Adelay.*(-1j*2*pi./cf.lambdak*cf.PU(2)/norm(cf.PU));

    
    cf.gain_mat = cf.gain*coe_mat; % size(alpha_TM)     NxK
    cf.dvec_mat = exp(-1j*2*pi/cf.lambdac*dis_vec.*cf.bse_coe);    % size(dvec_TM)   NxK
    Hs = cf.gain_mat.*dvec_mat.*Adelay;
    % cp.dvec_mat
    % cp.Adelay
    
%     D_H_x = Hs.*(-1j*2*pi./cf.lambdak).*((cf.PU(1)-cf.D_Rx_Loc(1,:))./vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1))' + ...
%          cf.rho*D_coe_x.*exp(-1j*2*pi*(vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1).')./cf.lambdak);
%     D_H_y = Hs.*(-1j*2*pi./cf.lambdak).*((cf.PU(2)-cf.D_Rx_Loc(2,:))./vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1))' + ...
%          cf.rho*D_coe_y.*exp(-1j*2*pi*(vecnorm(cf.D_Rx_Loc-cf.PU, 2, 1).')./cf.lambdak);
    D_H_x = cf.gain*D_coe_x.*(dvec_mat.*Adelay) + cf.gain_mat.*(D_dvec_x.*Adelay + dvec_mat.*D_Adelay_x);
    D_H_y = cf.gain*D_coe_y.*(dvec_mat.*Adelay) + cf.gain_mat.*(D_dvec_y.*Adelay + dvec_mat.*D_Adelay_y);
    
    D_H_rho = Hs/cf.rho;
    D_H_xi = 1j*Hs;


    for g = 1:cf.G
        Du_x = cf.WRx(:, g).'*D_H_x*cf.WTx(:, g).*(cf.X);
        Du_y = cf.WRx(:, g).'*D_H_y*cf.WTx(:, g).*(cf.X);
        Du_rho = cf.WRx(:, g).'*D_H_rho*cf.WTx(:, g).*(cf.X); 
        Du_xi = cf.WRx(:, g).'*D_H_xi*cf.WTx(:, g).*(cf.X);
        Du = [Du_x(:), Du_y(:), Du_rho(:), Du_xi(:)];
        FIM_S = FIM_S + 2/cf.sigma^2*real(Du'*Du);
    end

%     p.tBU_loc = get_dir_from_angle([cf.AOA 0]);
%     p.dBU = cf.d;
%     p.RotB = eye(3);
%     D_PU_theta = (p.tBU_loc(1)*p.RotB(:,2) - p.tBU_loc(2)*p.RotB(:,1))/(p.tBU_loc(1)^2 + p.tBU_loc(2)^2)/p.dBU;
%     D_PU_d = [cf.PU0; 0]/p.dBU;

    % FIM for state parameters
%     CRLB_S = (FIM_S)^-1;
    CRLB_S = get_EFIM_from_FIM(FIM_S,2)^-1;
    cf.PEB = sqrt(trace(CRLB_S(1:2, 1:2)));
    
    % FIM for channel parameters
    D_P_theta = cf.d*[-sind(cf.AOA) cosd(cf.AOA)]';
    D_P_d = [cosd(cf.AOA) sind(cf.AOA)]';
    J = eye(4);
    J(1:2, 1:2) = [D_P_theta(1:2) D_P_d(1:2)]';
    FIM_M = J*FIM_S*J.';
%     CRLB_M = (FIM_M,2)^-1;
    CRLB_M = get_EFIM_from_FIM(FIM_M,2)^-1;
    cf.AEB = rad2deg(sqrt(CRLB_M(1,1)));
    cf.DEB = sqrt(CRLB_M(2,2));
    

end