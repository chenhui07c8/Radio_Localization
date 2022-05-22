function output = alg_redundant_calculation(N, M)
    if M < 2 || M == N
        output = 0;
    else
        % case 1: N is even, M is even
        % case 2: N is even, M is odd. 0.
        % case 3: N is odd, M is even
        % case 4: N is odd, M is odd
        % redundant = (all - sym)/2
        if(mod(N,2)==0 && mod(M,2)==1) % no repeated
            output = (nchoosek(N, M) - 0)/2 + alg_redundant_calculation(N-1, M);
        else
            output = (nchoosek(N, M) - nchoosek(floor(N/2), floor(M/2)))/2 + alg_redundant_calculation(N-1, M);
        end
%         output = alg_redundant_calculation(N-1, M) + alg_redundant_calculation(N-1, M-1);
    end
    
end