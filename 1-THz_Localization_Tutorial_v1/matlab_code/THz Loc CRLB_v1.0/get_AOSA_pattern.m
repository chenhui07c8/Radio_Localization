function pattern = get_AOSA_pattern(beamforming_vec, temp_angle, c, params)
    TB = size(beamforming_vec, 2);

    pattern = zeros(TB, c.K);
    if params == "B"
        AE0 = c.Bae0;
        Nae = c.NaeB;
    elseif params == "M"
        AE0 = c.Mae0;
        Nae = c.NaeM;
    end
    % temp_angle = [0 0]';
    t_loc = get_dir_from_angle(temp_angle);  % get the BF direction vector for the b-th SA.
    for t_b = 1:TB
        coe = c.lambdak/c.lambdac;  % coefficient considering beam split (set to 1 to ignore beam split)
        tBFtemp = get_dir_from_angle(beamforming_vec(:, t_b));  % get the BF direction vector for the b-th SA.
        OmegaPsi = AE0'*t_loc - (AE0'*tBFtemp)*coe;
        AeqMat = 1/sqrt(Nae)*exp(1j*2*pi./c.lambdak.*(OmegaPsi));
        pattern(t_b, :) = abs(sum(AeqMat,1));  % imaginary part is 0 after summation
    end
    pattern = pattern/norm(pattern);
    % figure;plot(pattern)
end