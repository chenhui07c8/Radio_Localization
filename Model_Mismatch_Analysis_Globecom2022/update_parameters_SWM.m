function cp = update_parameters_SWM(cp)

    % HWI variables
    % cp.S = get_Tx_symbol_HW(cp);
    % PN (phase noise and carrier frequency offset)
    % size: K x G
    if cp.HW_PN == "True"
        cp.vec_PN_Tx = exp(1j*randn(cp.K, cp.G)*cp.sigma_PN_Tx);
        cp.vec_PN_Rx = exp(1j*randn(cp.K, cp.G)*cp.sigma_PN_Rx);
    elseif cp.HW_PN == "Tx"
        cp.vec_PN_Tx = exp(1j*randn(cp.K, cp.G)*cp.sigma_PN_Tx);
    elseif cp.HW_PN == "Rx"
        cp.vec_PN_Rx = exp(1j*randn(cp.K, cp.G)*cp.sigma_PN_Rx);
    else
        cp.vec_PN_Tx = ones(cp.K, cp.G);
    	cp.vec_PN_Rx = ones(cp.K, cp.G);
    end
    
    % CFO (carrier frequency offset)
    % size: K x G
    if cp.HW_CFO == "True"
        cp.vec_CFO_Tx = exp(1j*2*pi*cp.epsilon_CFO_Tx/cp.N_Tx*( (0:(cp.K-1))' - (cp.K-1)/2 ));
        cp.vec_CFO_Rx = exp(1j*2*pi*cp.epsilon_CFO_Rx/cp.N_Rx*( (0:(cp.K-1))' - (cp.K-1)/2 ));
    elseif cp.HW_CFO == "Tx"
        cp.vec_CFO_Tx = exp(1j*2*pi*cp.epsilon_CFO_Tx/cp.N_Tx*( (0:(cp.K-1))' - (cp.K-1)/2 ));
        cp.vec_CFO_Rx = ones(cp.K, 1);
    elseif cp.HW_CFO == "Rx"
        cp.vec_CFO_Tx = ones(cp.K, 1);
        cp.vec_CFO_Rx = exp(1j*2*pi*cp.epsilon_CFO_Rx/cp.N_Rx*( (0:(cp.K-1))' - (cp.K-1)/2 ));   
    else
        cp.vec_CFO_Tx = ones(cp.K, 1);
        cp.vec_CFO_Rx = ones(cp.K, 1);
    end
    
    % MC
    C_Tx = toeplitz([1, cp.u_MC_Tx, zeros(1, (cp.N_Tx)-length(cp.u_MC_Tx)-1)]);
    C_Rx = toeplitz([1, cp.u_MC_Rx, zeros(1, (cp.N_Rx)-length(cp.u_MC_Rx)-1)]);

    if cp.HW_MC == "True"
        cp.C_Tx = C_Tx(1:cp.N_Tx, 1:cp.N_Tx) + (randn(cp.N_Tx) + 1j*randn(size(cp.N_Tx)))*(cp.sigma_MC_Tx/sqrt(2));
        cp.C_Rx = C_Rx(1:cp.N_Rx, 1:cp.N_Rx) + (randn(cp.N_Rx) + 1j*randn(size(cp.N_Rx)))*(cp.sigma_MC_Rx/sqrt(2));
    elseif  cp.HW_MC == "Tx"
        cp.C_Tx = C_Tx(1:cp.N_Tx, 1:cp.N_Tx) + (randn(cp.N_Tx) + 1j*randn(size(cp.N_Tx)))*(cp.sigma_MC_Tx/sqrt(2));
        cp.C_Rx = C_Rx(1:cp.N_Rx, 1:cp.N_Rx);
    elseif  cp.HW_MC == "Rx"
        cp.C_Tx = C_Tx(1:cp.N_Tx, 1:cp.N_Tx);
        cp.C_Rx = C_Rx(1:cp.N_Rx, 1:cp.N_Rx) + (randn(cp.N_Rx) + 1j*randn(size(cp.N_Rx)))*(cp.sigma_MC_Rx/sqrt(2));
    else
%         cp.C_Tx = eye((cp.N_Tx));
%         cp.C_Rx = eye((cp.N_Rx));
        cp.C_Tx = C_Tx(1:cp.N_Tx, 1:cp.N_Tx);
        cp.C_Rx = C_Rx(1:cp.N_Rx, 1:cp.N_Rx);
    end
    
    % AR: array response calibration
    % size: N x L
    if cp.HW_AR == "True"
        cp.vec_B_Tx = ones(cp.N_Tx, cp.L) + (randn(cp.N_Tx, cp.L) + 1j*randn(cp.N_Tx, cp.L))*(cp.sigma_AR_Tx/sqrt(2));
        cp.vec_B_Rx = ones(cp.N_Rx, cp.L) + (randn(cp.N_Rx, cp.L) + 1j*randn(cp.N_Rx, cp.L))*(cp.sigma_AR_Rx/sqrt(2));
    elseif cp.HW_AR == "Tx"
        cp.vec_B_Tx = ones(cp.N_Tx, cp.L) + (randn(cp.N_Tx, cp.L) + 1j*randn(cp.N_Tx, cp.L))*(cp.sigma_AR_Tx/sqrt(2));
    elseif cp.HW_AR == "Rx"
        cp.vec_B_Rx = ones(cp.N_Rx, cp.L) + (randn(cp.N_Rx, cp.L) + 1j*randn(cp.N_Rx, cp.L))*(cp.sigma_AR_Rx/sqrt(2));
    else
        cp.vec_B_Tx = ones(cp.N_Tx, cp.L);
        cp.vec_B_Rx = ones(cp.N_Rx, cp.L);
    end
    
    % AG: array geometry calibration
    if cp.HW_AG == "True"
        delta_D_Tx = randn(size(cp.D_Tx))*cp.sigma_AG_Tx;
        delta_D_Rx = randn(size(cp.D_Rx))*cp.sigma_AG_Rx;
    elseif cp.HW_AG == "Tx"
        delta_D_Tx = randn(size(cp.D_Tx))*cp.sigma_AG_Tx;
        delta_D_Rx = zeros(size(cp.D_Rx));
    elseif cp.HW_AG == "Rx"
        delta_D_Tx = zeros(size(cp.D_Tx));
        delta_D_Rx = randn(size(cp.D_Rx))*cp.sigma_AG_Rx;
    else
        delta_D_Tx = zeros(size(cp.D_Tx));
        delta_D_Rx = zeros(size(cp.D_Rx));
    end
    cp.D_Tx = cp.D_Tx + delta_D_Tx;
    cp.D_Rx = cp.D_Rx + delta_D_Rx;
    
    % IQI
    if cp.HW_IQI == "True"
        cp.IQI_alpha_Tx = 0.5*(1+cp.IQI_m_Tx*exp(1j*cp.IQI_phi_Tx));
        cp.IQI_beta_Tx = 0.5*(1-cp.IQI_m_Tx*exp(1j*cp.IQI_phi_Tx));
        cp.IQI_alpha_Rx = 0.5*(1+cp.IQI_m_Rx*exp(1j*cp.IQI_phi_Rx));
        cp.IQI_beta_Rx = 0.5*(1-cp.IQI_m_Rx*exp(1j*cp.IQI_phi_Rx));
    elseif cp.HW_IQI == "Tx"
        cp.IQI_alpha_Tx = 0.5*(1+cp.IQI_m_Tx*exp(1j*cp.IQI_phi_Tx));
        cp.IQI_beta_Tx = 0.5*(1-cp.IQI_m_Tx*exp(1j*cp.IQI_phi_Tx));
    elseif cp.HW_IQI == "Rx"
        cp.IQI_alpha_Rx = 0.5*(1+cp.IQI_m_Rx*exp(1j*cp.IQI_phi_Rx));
        cp.IQI_beta_Rx = 0.5*(1-cp.IQI_m_Rx*exp(1j*cp.IQI_phi_Rx));
    else
    	cp.IQI_alpha_Tx = 1;
        cp.IQI_beta_Tx = 0;
        cp.IQI_alpha_Rx = 1;
        cp.IQI_beta_Rx = 0;
    end


    % 3: received symbol
    cp = get_Rx_symbol_HWI(cp);

end