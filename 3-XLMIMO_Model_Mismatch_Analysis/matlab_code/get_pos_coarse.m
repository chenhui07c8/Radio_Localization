function initial = get_pos_coarse(y, cf)

% coarse: angle
    angle_all = cf.AOA_candidate;
    cost_ml = zeros(1,length(angle_all));
    for ang_i = 1:length(angle_all)
        theta = angle_all(ang_i);
        cf.AOA = theta;
%         c2.AOA = 59.0362;
        cost = get_cost_angle_coarse(y, cf);
        cost_ml(ang_i) = cost;
    end
    % figure;plot(angle_all, cost_ml);
    [a,b] = min((cost_ml));
    coarse_aoa = angle_all(b);

% coarse: delay
    d_all = cf.d_candidate;
    cost_ml = zeros(1,length(d_all));
    for d_i = 1:length(d_all)
        cf.AOA = coarse_aoa;
%         c2.AOA = 59.0362;
        cf.d = d_all(d_i);
        cost = get_cost_delay_coarse(y, cf);
        cost_ml(d_i) = cost;
    end
    % figure;plot(d_all, cost_ml);
    [a,b] = min((cost_ml));
    coarse_d = d_all(b);

    
    initial = coarse_d*get_dir_from_angle([coarse_aoa 0]');

        
end