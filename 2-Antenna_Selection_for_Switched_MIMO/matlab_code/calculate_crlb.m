function crlb = calculate_crlb(d, snrdB)
    result = 0;
    M = length(d);
    snr = M*10.^(snrdB/10);
    for i = 1:M
        result = result + (d(i)-mean(d))^2;
    end
    result = result/M;

    CRB = 1/2/pi^2/result./snr;
    % crlb = CRB*(180/pi/cos(39.5/180*pi))^2;
    crlb = CRB;


end