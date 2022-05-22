% Author: hui.chen@kaust.edu.sa
% last update: 24-May-2021
% reference: 
% difference in weights for different directions.

% M: m-by-3 BS positions
% A: m-by-3 direction vectors
% Err_Dir;  % STD of the DOA (sigma)

function loc_est = getLineLineIntersection(M, A, DirVariance)

    [m, n] = size(M);  % m base stations
    C = nchoosek(1:m,2);
    [m_c, n_c] = size(C);   % m_c, number of array pairs
    sjk = zeros(m_c, 3);
    skj = zeros(m_c, 3);


    for i = 1:m_c
        j = C(i,1);
        k = C(i,2);
        mj = M(j,:);
        mk = M(k,:);
        aj = A(j,:);
        ak = A(k,:);    
        djk = abs(dot(cross(aj,ak),(mj-mk)))/norm(cross(aj,ak));
        b = mk-mj-djk*(cross(aj,ak));   % overconstrained matrix equation
        Amat = [aj' -ak'];
        rjk = (Amat'*Amat)^-1*Amat'*b';   % r is estimated distance using pure angles
        rj = rjk(1);
        rk = rjk(2);
        sjk(i,:) = mj + rj*aj;
        skj(i,:) = mk + rk*ak;
    end


    % weighting....of linear interpolation estimator
    W = zeros(1,m_c);
    S = zeros(3,m_c);

    for pair = 1:m_c
        W(pair) = 1/DirVariance(C(pair, 1))/DirVariance(C(pair, 2));  % weight of each pair
        S(:, pair) = (1/DirVariance(C(pair, 1))*sjk(pair,:)' + 1/DirVariance(C(pair, 2))*skj(pair,:)')...
            /(1/DirVariance(C(pair, 1)) + 1/DirVariance(C(pair, 2)));  % weighted position of each pair
    end


    % weighted results
    loc_est = sum(W.*S, 2);


end