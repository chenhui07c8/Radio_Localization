function output = alg_unique_calculation(N, M)
    if M == N
        output = 1;
    elseif M == 0
        output = 0;
    else
        output = alg_unique_calculation(N-1, M) + alg_unique_calculation(N-1, M-1);
    end
    
end