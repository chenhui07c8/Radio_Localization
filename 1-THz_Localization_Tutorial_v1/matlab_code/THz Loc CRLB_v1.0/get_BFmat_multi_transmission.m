function cMul = get_BFmat_multi_transmission(c)
    
    T_B = c.T_B;
    T_MR = c.T_MR;
    T_MB = c.T_MB;
    angle_max = 80;
    beams_B = [linspace(-angle_max, angle_max, T_B); zeros(1, T_B)];
    beams_MB = [linspace(-angle_max, angle_max, T_MB); zeros(1, T_MB)];
    beams_MR = [linspace(-angle_max, angle_max, T_MR); zeros(1, T_MR)];

    cMul.T = T_MB*T_B + T_MR;    % T_B rounds for UE, 1 round for RIS. (1 round = T_M) 

    % activation matrix (Antenna Array)
    cMul.NR_dim_mat = repmat(c.NR_dim, [1, cMul.T]);
    cMul.NR_dim_mat(:, 1:T_MB*T_B) = 0;
    cMul.NB_dim_mat = repmat(c.NB_dim, [1, cMul.T]);
    cMul.NM_dim_mat = repmat(c.NM_dim, [1, cMul.T]);

    % beamforming matrix (Beamforming Angle)
    cMul.BFmatB = zeros(2, c.NB, cMul.T);
    cMul.BFmatM = zeros(2, c.NM, cMul.T);

    % BM channel: t_b = 1:T_B; BRM channel: t_b = T_B + 1;
    for t_b = 1:T_B
        for t_mb = 1:T_MB
            cMul.BFmatB(:, :, T_MB*(t_b-1) + t_mb) = repmat(beams_B(:, t_b), [1, c.NB]);
            cMul.BFmatM(:, :, T_MB*(t_b-1) + t_mb) = repmat(beams_MB(:, t_mb), [1, c.NM]);
        end
    end
    % BRM channel
    for t_mr = 1:T_MR  
        cMul.BFmatB(:, :, T_B*T_MB + t_mr) = repmat([c.phiBR_loc 0]', [1, c.NB]);
        cMul.BFmatM(:, :, T_B*T_MB + t_mr) = repmat(beams_MR(:, t_mr), [1, c.NM]);
    end
    
end