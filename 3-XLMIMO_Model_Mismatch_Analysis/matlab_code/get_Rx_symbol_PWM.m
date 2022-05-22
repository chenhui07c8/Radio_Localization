function c = get_Rx_symbol(c)
    % Transmitted symbols before transmission
    c.AstTx = exp(1j*pi*c.D_Tx*sind(c.theta_AOD));  % steering vector at the Rx
    c.AstRx = exp(1j*pi*c.D_Rx*sind(c.theta));  % steering vector at the Rx
    c.H = c.rho*c.AstRx*c.AstTx.';
    c.Adelay = exp(-1j*2*pi*(0:c.K-1)'*c.delta_f*c.d0/c.c);
    % c.Adoppler = exp(1j*2*pi*(1:c.G-1)*c.Tsym*c.v0/c.lambdac);

    c.u = zeros(c.K, c.G);
    for g = 1:c.G
        c.u(:, g) = c.WRx(:, g).'*c.H*(c.Adelay.*c.s).'; % c.G*c.K
    end

end