function cp = get_Rx_symbol(cp)
%     cp.AOD = atan2d(cp.PU(2), cp.PU(1));    % angle of departure
%     cp.d = norm(cp.PU);

    
    u0 = zeros(cp.M_Rx, cp.K, cp.G);
    
	% get channels
    dis_Rx = -cp.lambdac/2*(cp.D_Rx-mean(cp.D_Rx))*sind(cp.AOA);
    cp.AstRx = exp(-1j*2*pi/cp.lambdac*dis_Rx);  % steering vector at the Rx
    cp.Adelay = exp(-1j*2*pi*cp.fk*cp.d/cp.c);
    cp.AstTx = 1;
%     exp(1j*pi*sind(cp.AOA)*(2*(1:cp.N)'-cp.N-1)/2)
    % AR
    if cp.wave_type == "PWM"
        % PWM
        H = cp.gain*cp.AstRx.*cp.Adelay;
        H1 = H;
%                 H = cp.gain*cp.AstRx*(cp.AstTx).'.*cp.Adelay;
    elseif cp.wave_type == "SWM"
        if cp.NF_SNS == "True"
            coe_mat = cp.lambdak./cp.lambdac.*(norm(cp.PU)./vecnorm(cp.D_Rx_Loc-cp.PU, 2, 1)).';
        else
            coe_mat = ones(cp.N, cp.K);
        end

        if cp.NF_SWM == "True"
            dis_vec = vecnorm(cp.D_Rx_Loc-cp.PU, 2, 1).' - norm(cp.PU);
        else
            dis_vec = dis_Rx;
        end

        if cp.NF_BSE == "True" 
            bse_coe = cp.lambdac./cp.lambdak;
        else
            bse_coe = ones(size(cp.fk));
        end
        cp.gain_mat = cp.gain*coe_mat; % size(alpha_TM)     NxK
        cp.dvec_mat = exp(-1j*2*pi/cp.lambdac*dis_vec.*bse_coe);    % size(dvec_TM)   NxK
        % cp.dvec_mat(:,1)
        H = cp.gain_mat.*cp.dvec_mat.*cp.Adelay;
        H2 = H;
        cp.bse_coe = bse_coe;

    end
    % get received symbols
    for g = 1:cp.G
            cp.X = sqrt(cp.P)*cp.s.';

%             norm(H1-H2, 'fro')
            % exp(1j*pi*cp.D*sind(cp.theta));
            % equals to: H = c.gain*diag(c.vec_B_Rx)*exp(1j*pi*c.D*sind(c.theta));
            if cp.array_structure == "Digital"
                u0(:, :, g) = H.*(cp.X);
            elseif cp.array_structure == "Analog"
                u0_k = cp.WRx(:, g).'*H.*(cp.WTx(:, g).*(cp.X));    % for different subcarriers
%                 u0_k = cp.WRx(:, g).'*H.*((cp.X));    % for different subcarriers

                u0(:, :, g) = u0_k; % c.G*c.K
            end
    end
    
    cp.u = u0;
end


