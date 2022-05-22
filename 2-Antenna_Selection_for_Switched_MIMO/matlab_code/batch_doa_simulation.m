% doa_center = simulated target DOA;
% thetam: searching range of MLE from -thetam to thetam
% S_lambda: sensor location in lambda/2 
% snr: in dB
% sim_times: simulation times.

function [errors, results, estR] = batch_doa_simulation(target_all, theta_max, dv, snr)
    Ns = 1;  % number of snapshot
    M = length(dv); % number of sensors
    errors = zeros(1, length(target_all)); % DOA estimation error
    results = zeros(1, length(target_all)); % DOA estimation results
    for i_angle = 1:length(target_all)
        temp_angle = target_all(i_angle);

        % calculate covariance matrix
        theta   = temp_angle/180*pi;       
        a      = exp(1j*pi*dv*sin(theta))';
        a      = sqrt(M)*a/norm(a); 
        s      = ones(Ns, 1);
        y0     = a*transpose(s);  
        y      = awgn(y0, snr, 'measured');
        estR   = zeros(M, M);
        for k=1:Ns
            estR = estR + y(:, k)*y(:, k)';
        end
        output = alg_mle(estR, dv, theta_max, [0.2 0.01]);
        results(i_angle) = output;
        errors(i_angle) = sin(output/180*pi) - sin(temp_angle/180*pi);
    end

end
