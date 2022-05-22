function output = alg_nn(target_doa, snr, variables, D, N)
    
    a = [target_doa, snr];
    for i = 1:length(variables)/2
        l = a*variables{2*i-1} + variables{2*i};
        a = max(0, l);
    end
%     output = 1 ./ (1 + exp(-a));
    [B, I] = sort(-a);
%     output = zeros(1,length(D));
%     output(I(1:N)) = 1;
%     output = I(1:N);
%     
    output = sort(D(I(1:N)));
end