function [FIM, FIMmat] = get_FIM_multi_transmission(c, cM)

    % [uBmat = FIMmat]
%     if c.pos_type + c.ori_type == "2D1D"
%         fim_size = 10;
%     elseif c.pos_type + c.ori_type == "3D3D"
%         fim_size = 14;
%     end

%     uBmat = zeros(c.NB, c.K, cM.T);
%     FIMmat = zeros(fim_size, fim_size, cM.T);

%     c.G = 1;
    
    [FIM0] = get_fim_uplink(c);
    FIMmat = zeros(size(FIM0,1), size(FIM0,2), c.G);
    c.symbol_type = "Customized"; % Ones, Random, RandomSymbol, Customized
    c.RIS_optimizer = "Customized";
    for t = 1:c.G
%         c.NB_dim = cM.NB_dim_mat(:,t);  
%         c.NR_dim = cM.NR_dim_mat(:,t);  
%         c.NM_dim = cM.NM_dim_mat(:,t);  
        c.rng_ind = t;
        c.BFmatB = cM.BFmatB(:,:, t);
        c.BFmatM = cM.BFmatM(:,:, t);
        c.Xmk = cM.XmkMat(:,:, t);
        if(c.LR>0)
            c.Omega = cM.OmegaMat(:,t);
            c.BFmatRB = cM.BFmatRB(:,:, t);
            c.BFmatRM = cM.BFmatRM(:,:, t);
        end       
        c = updateParameters(c);
        [FIM0] = get_fim_uplink(c);
        FIMmat(:,:,t) = FIM0;
    end
    FIM = sum(FIMmat, [3]);
%     FIM = cell2mat(FIMmat);

    FIM = FIM/c.G;
end