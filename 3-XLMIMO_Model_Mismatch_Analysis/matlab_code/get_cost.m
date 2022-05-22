function cost = get_cost(y, cf)
    if cf.wave_type == "PWM"
        u = cf.u/cf.gain;
        y0 = y(:);
        u0 = u(:);
        gain = (u0)'*(y0)/norm(u0)^2;
%         cost = norm(y0 - gain*u0);

    elseif cf.wave_type == "SWM"
%         coe_mat = cp.fk./cp.fc.*(norm(cp.PU)./vecnorm(cp.D_Rx_Loc-cp.PU, 2, 1)).';
%         u = cf.u./coe_mat;
        u = cf.u/cf.gain;
        y0 = y(:);
        u0 = u(:);
        gain = (u0)'*(y0)/norm(u0)^2;
%         cost = norm(y0 - gain*u0);

    end
    cost = norm(y0 - gain*u0);
end