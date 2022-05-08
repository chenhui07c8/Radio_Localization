function cost = get_cost_angle_coarse(y, cf)

    cf = get_Rx_symbol(cf);
    u = cf.u./(cf.gain*cf.Adelay);
    % mean(diff( phase(cf.u(:,1))))
    rho_vec = zeros(1, cf.K);
    for ki = 1:cf.K
        y0 = y(:, ki);
        u0 = u(:, ki);
        rho_vec(ki) = (u0)'*(y0)/norm(u0)^2;
    end
    hat_y = rho_vec.*u;
    cost = norm(y(:) - hat_y(:));

%     cost = norm(y - rho_vec.*u, 'fro');
    %           cp.Adelay = exp(-1j*2*pi*(0:cp.K-1)'*cp.delta_f*cp.d/cp.c);
    %           mean(diff(phase(rho_vec)))/(2*pi*cp.delta_f/cp.c)
    %           mean(diff(phase(c2.Adelay)))

end