function [S, x0] = linear_amb(D, M, delta)

N = length(D);
K = M;
n = -(N-1)/2:1:(N-1)/2;
n = n/2;
id = 1:N;
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;

%define the array position matrix
dx = d*n';
dxx = dx.*dx;

%avoiding ambiguity
vjs1 = exp(1i*k0*d*n'*0.12);
vjs2 = exp(1i*k0*d*n'*0.24);
vjs3 = exp(1i*k0*d*n'*0.36);
vjs4 = exp(1i*k0*d*n'*0.48);
W = [vjs1,vjs2,vjs3,vjs4];

id = 1:N;
V = nchoosek(id,K);
num = nchoosek(N,K);

e = ones(N,1);
mu = 10;
Q = 4;
x0 = [ones(K/2,1);zeros(N-K,1);ones(K/2,1)];
for q = 1:Q
    cvx_begin quiet
        variable x(N);
        minimize (-x'*dxx+mu*(e'*x-2*x0'*x));
        subject to
        0 <= x'*dx <= 0;
        K <= e'*x <= K;
        0 <= x <= 1;
        for i = 1:size(W,2)
            W = (1/(K*K))*real(W(:,i)*W(:,i)');
            x'*W*x <= delta;
        end
    cvx_end
    x0 = x;
end
[B, I] = sort(x0, 'descend');
S = D(sort(I(1:K)));

end