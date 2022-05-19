function [output, index, value] = alg_threshold_greedy(target_doa, ref_doa, snrdB, D, M)
    d = D;
    L_start = length(D);
    index = ones(size(D));
    if(length(target_doa) == 1)    % single point
        [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, D', snrdB);
        optimal = Es(2);
        for loop_d = 1:(L_start - M)
    %         loop_d = 2
            dcan = zeros(1, length(d));
            for i = 1:length(d)
                d_temp = [d(1:(i-1)) d((i+1):end)];
                [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, d_temp', snrdB);
                dcan(i) = Es(2);
            end

            [B, I] = sort(dcan);
            i = I(1);
            index(D==d(i)) = 0;
            d = [d(1:(i-1)) d((i+1):end)];
        end
        value = B(1) - optimal;
        % last sensor value
        output = d;
    else    % averaging multiple DOAs
        optimal_arr = zeros(1,length(target_doa));
        for doa_i = 1:length(target_doa)
            [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa(doa_i), ref_doa, D', snrdB);
            optimal_arr(doa_i) = Es(2);
        end
        optimal = min(optimal_arr);
        
        for loop_d = 1:(L_start - M)
    %         loop_d = 2
            dcan = zeros(1, length(d));
            for i = 1:length(d)
                d_temp = [d(1:(i-1)) d((i+1):end)];
                dcan_arr = zeros(1, length(target_doa));
                for doa_i = 1:length(target_doa)
                    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa(doa_i), ref_doa, d_temp', snrdB);
                    dcan_arr(doa_i) = Es(2);
                end
                dcan(i) = max(dcan_arr);    % worst case
                % dcan(i) = mean(dcan_arr); % average
                % dcan(i) = min(dcan_arr);  % best case

            end

            [B, I] = sort(dcan);
            i = I(1);
            index(D==d(i)) = 0;
            d = [d(1:(i-1)) d((i+1):end)];
        end
        value = B(1) - optimal;
        % last sensor value
        output = d;
    end
        
        
end