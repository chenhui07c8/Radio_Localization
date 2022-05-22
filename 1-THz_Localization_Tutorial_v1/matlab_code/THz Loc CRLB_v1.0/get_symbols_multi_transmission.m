function uBmat = get_symbols_multi_transmission(c, cM)

    % [uBmat = FIMmat]
%     if c.pos_type + c.ori_type == "2D1D"
%         fim_size = 10;
%     elseif c.pos_type + c.ori_type == "3D3D"
%         fim_size = 14;
%     end

    uBmat = zeros(c.NB, c.K, c.G);
%     FIMmat = zeros(fim_size, fim_size, cM.T);

    for t = 1:c.G
%         c.NB_dim = cM.NB_dim_mat(:,t);  
%         c.NR_dim = cM.NR_dim_mat(:,t);  
%         c.NM_dim = cM.NM_dim_mat(:,t);  

        c.BFmatB = cM.BFmatB(:,:, t);
        c.BFmatM = cM.BFmatM(:,:, t);
        c = updateParameters(c, "All");
        c.Xmk = cM.XmkMat(:,:, t);
        c.Omega = get_ris_phase(c);   % get ris phase (w/o beamsplit, c.beamsplit)
        c = get_Rx_symbols(c);
        uBmat(:, :, t) = c.uB;
%         [FIM0] = get_fim(c);
%         FIMmat(:, :, t) = FIM0;
    end
%     FIM = sum(FIMmat, [3]);


end