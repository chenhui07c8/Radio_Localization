function [S, x0] = linear_iso_search(D,M,delta)


    N = length(D);
    K = M;
    n = -(N-1)/2:1:(N-1)/2;
    n = n/2;
    id = 1:N;
    lambda = 1;
    d = lambda/2;
    k0 = 2*pi/lambda;

    %define the array position matrix (2N multiplications)
    dx = d*n';
    dxx = dx.*dx;

    %avoiding ambiguity (4N multiplications)
    vjs1 = exp(1i*k0*d*n'*0.12);
    vjs2 = exp(1i*k0*d*n'*0.24);
    vjs3 = exp(1i*k0*d*n'*0.36);
    vjs4 = exp(1i*k0*d*n'*0.48);
    W = [vjs1,vjs2,vjs3,vjs4];

    id = 1:N;
    V = nchoosek(id,K);
    num = nchoosek(N,K);    
    e = ones(K,1);
    %searching
    y = zeros(num,1);
    for i = 1:num   % N choose K iterations: F_NM*M + \F_NM*4M
        if (e'*dx(V(i,:))==0)
            flag = 0;
            for k = 1:size(W,2)
                v = W(:,k);
                if ((1/K)*abs(v(V(i,:))'*e) <= delta)
                    flag = flag + 1;
                end
            end
            if (flag == size(W,2))
                y(i) = e'*dxx(V(i,:));
            end
        end  
    end
    [ymax,Imax] = max(y);
    x0 = zeros(N,1);
    for i = 1:K
    x0(V(Imax,i)) = 1;
    end
    [B, I] = sort(x0, 'descend');
    S = D(sort(I(1:K)));
end



