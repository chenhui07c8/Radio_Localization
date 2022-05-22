function output = alg_crlb_greedy(snrdB, D, N)
    d = D;
    L_start = length(D);
    for loop_d = 1:(L_start - N)
        dcan = zeros(1, length(d));
        for i = 1:length(d)
            d_temp = [d(1:(i-1)) d((i+1):end)];
%             [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, d_temp', snrdB);
            crlb = calculate_crlb(d_temp, snrdB);
            dcan(i) = crlb;
        end
        [B, I] = sort(dcan);
        i = I(1);
        d = [d(1:(i-1)) d((i+1):end)];
%         d;
    end

    output = d;
end