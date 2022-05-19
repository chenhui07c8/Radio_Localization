function output = alg_mle(estR, S_lambda, thetam, resolution)
% refine into 0.0001;
    Nr = length(S_lambda);
    angle_all = -thetam:resolution(1):thetam;
    cost_ml = zeros(1,length(angle_all));
    for ang_i = 1:length(angle_all)
        theta = angle_all(ang_i);
        av = exp(1j*pi*S_lambda*sin(theta/180*pi))';
%         A = 1-av*(av'*av)^-1*av';
%         cost_ml(ang_i) = trace(A*estR);
        cost_ml(ang_i) = av'*estR*av/Nr;
    end
    [~,b] = max(abs(cost_ml));
    coarse = angle_all(b);
    % figure;plot(angle_all, abs(cost_ml))
    
    angle_all = (-resolution(1):resolution(2):resolution(1)) + coarse;
    cost_ml = zeros(1,length(angle_all));
    for ang_i = 1:length(angle_all)
        theta = angle_all(ang_i);
        av = exp(1j*pi*S_lambda*sin(theta/180*pi))';
%         A = 1-av*(av'*av)^-1*av';
%         cost_ml(ang_i) = trace(A*estR);
        cost_ml(ang_i) = av'*estR*av/Nr;      
    end
    [a,b] = max(abs(cost_ml));
    fine = angle_all(b);
    output = fine;
%     angle_all = (-0.02:0.001:0.02) + fine;
%     costml = zeros(1,length(angle_all));
%     for ang_i = 1:length(angle_all)
%         ttt = angle_all(ang_i);
%         av = exp(1j*pi*S_lambda*sin(ttt/180*pi))';
%         av = sqrt(Nr)*av/norm(av);       
%         costml(ang_i) = av'*estR*av;
%     end
%     [a,b] = max(abs(costml));
%     output = angle_all(b);


end




