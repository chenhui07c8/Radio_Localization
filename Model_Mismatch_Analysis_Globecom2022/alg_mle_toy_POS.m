% function: get the AOD estimation
% u0: observed data
% c: channel information
% resolution: searching space (coarse fine)

function [PU_est, gain_optimal] = alg_mle_toy_POS(y, cf)


    iter_max = 100;
    PU_iter = zeros(2, iter_max);
    cost_iter = zeros(1, iter_max);

%     pos_optimal = [4 4.3]';
%     pos_optimal = c2.PU0 + 0.5*[1 1]';
    pos_optimal = cf.PU;

    gain_optimal = 1;
    cost_optimal = 1e9;
%         results_cell{i}(2,j) = (theta - doa);
    step_size = 0.2;
%     step_dir = [-1 -1]'/sqrt(2);
    step_dir = [0 0]';
%         grad_x = 0.3;
%         grad_y = 0.3;
%         c2 = cf;
    for iter_i = 1:iter_max
        % next step; iter_i = 2
        cf.PU = pos_optimal - step_size*step_dir;
        cf.AOD = 0;    % angle of arrival
        cf.AOA = atan2d(cf.PU(2), cf.PU(1));    % angle of departure
        cf.d = norm(cf.PU);
        cf = get_Rx_symbol(cf);
        u = cf.u/cf.gain;
        y0 = y(:);
        u0 = u(:);
        gain = (u0)'*(y0)/norm(u0)^2;
        cost_temp1 = norm(y0 - gain*u0);

        % check cost
        
        
        % check next cost
        if(cost_optimal <= cost_temp1)
            step_size = step_size/2;
        else
            pos_optimal = pos_optimal - step_size*step_dir;
            gain_optimal = gain;
            cost_optimal = cost_temp1;
            % x axis: numerical gradient (could be changed with closed-form)
            cf.PU = pos_optimal + [1e-7; 0];
            cf.AOD = 0;    % angle of arrival
            cf.AOA = atan2d(cf.PU(2), cf.PU(1));    % angle of departure
            cf.d = norm(cf.PU);
            cf = get_Rx_symbol(cf);
            u = cf.u/cf.gain;
            y0 = y(:);
            u0 = u(:);
            gain = (u0)'*(y0)/norm(u0)^2;
            cost_temp_x = norm(y0 - gain*u0);                
            step_dir(1) = (cost_temp_x - cost_optimal);
            % y axis
            cf.PU = pos_optimal + [0; 1e-7];
            cf.AOD = 0;    % angle of arrival
            cf.AOA = atan2d(cf.PU(2), cf.PU(1));    % angle of departure
            cf.d = norm(cf.PU);
            cf = get_Rx_symbol(cf);
            u = cf.u/cf.gain;
            y0 = y(:);
            u0 = u(:);
            gain = (u0)'*(y0)/norm(u0)^2;
            cost_temp_y = norm(y0 - gain*u0);                
            step_dir(2) = (cost_temp_y - cost_optimal);
            step_dir = step_dir/norm(step_dir);
        end
        PU_iter(:, iter_i) = pos_optimal;
        cost_iter(iter_i) = cost_optimal;
    end
    PU_est = PU_iter(:, end);
% figure;plot(cost_iter)

end




