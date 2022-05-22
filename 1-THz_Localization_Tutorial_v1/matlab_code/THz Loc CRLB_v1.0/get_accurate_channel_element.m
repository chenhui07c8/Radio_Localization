function output = get_accurate_channel_element(c, b, m)

    H_sa = zeros(c.NsaB, c.NsaM, length(c.fk)); % channel matrix between SA-b and SA-m
    Bae = c.PB + c.RotB*(c.Bae0 + c.B0(:,b));
    Mae = c.PM + c.RotM*(c.Mae0 + c.M0(:,m));

    for bar_b = 1:c.NsaB
        for bar_m = 1:c.NsaM
            dbm = norm(Bae(:, bar_b)-Mae(:, bar_m));    % delay to the Rx
            rho = c.lambdac/4/pi/dbm;
            xi = exp(-1j*2*pi/c.lambdac*dbm);
            H_sa(bar_b, bar_m, :) = rho.*xi;
        end
    end

    % analog beamforming (at SA)
    tb_loc = get_dir_from_angle(c.BFmatB(:, b));
    tm_loc = get_dir_from_angle(c.BFmatM(:, m));
    
    beamformer_BS = exp(-1j*2*pi/c.lambdac*transpose(c.Bae0)*tb_loc);
    beamformer_MD = exp(-1j*2*pi/c.lambdac*transpose(c.Mae0)*tm_loc);

    channel_elements = zeros(1, length(c.fk));
    for i = 1:length(c.fk)
        H_temp = H_sa(:,:,i);
        channel_elements(i) = 1/sqrt(c.NsaB*c.NsaM)*transpose(beamformer_BS)*H_temp*beamformer_MD;
    end
    output = channel_elements;

end