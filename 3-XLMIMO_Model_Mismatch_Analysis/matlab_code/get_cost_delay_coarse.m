function  cost = get_cost_delay_coarse(y, cf)


%             cf = get_Rx_symbol_HWI(cf);
%             udelay = cf.Adelay;
    udelay = exp(-1j*2*pi*(0:cf.K-1)'*cf.delta_f*cf.d/cf.c);
    rho_vec = zeros(cf.M_Rx, cf.G);
%     for gi = 1:cf.G
%         y0 = y(:, :, gi).';
%         u0 = udelay;
%         rho_vec(gi) = (u0)'*(y0)/norm(u0)^2;
%     end
    for mi = 1:cf.M_Rx
        y0 = y(mi, :).';
        u0 = udelay;
        rho_vec(mi) = (u0)'*(y0)/norm(u0)^2;
    end
    hat_y = rho_vec.*udelay.';
    cost = norm(y(:) - hat_y(:));

%     cost = norm(y - rho_vec.*udelay.', 'fro');

end

