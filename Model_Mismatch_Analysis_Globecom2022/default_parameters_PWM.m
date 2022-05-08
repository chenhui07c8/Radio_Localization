function cp = default_parameters_HWI(cp)
    % PN: diag(exp(1j*randn(c.N, 1)*sigma_PN))
    cp.sigma_PN_Tx = 0.02;
    cp.sigma_PN_Rx = 0.02;
    
    % CFO: diag(exp(1j*2*pi*epsilon/c.N*( (0:(N-1)) - (N-1)/2) ) + residue
%     cp.epsilon_CFO_Tx = 0;
%     cp.epsilon_CFO_Rx = 0;
    cp.sigma_CFO_Tx = 1e-2;
    cp.sigma_CFO_Rx = 1e-2;
    cp.epsilon_CFO_Tx = randn(1,1)*cp.sigma_CFO_Tx;
    cp.epsilon_CFO_Rx = randn(1,1)*cp.sigma_CFO_Rx;
    % debug value
%     cp.epsilon_CFO_Tx = 1;
%     cp.epsilon_CFO_Rx = 1;
%     cp.sigma_CFO_Tx = 0;
%     cp.sigma_CFO_Rx = 0;
    
    
    % MC: toeplitz(cp.u_MC_Tx) + (randn(c.N)+1j*randn(c.N))*sigma_MC/sqrt(2)
    cp.sigma_MC_Tx = 0.01;
    cp.sigma_MC_Rx = 0.01;

    % debug value
%     cp.u_MC_Tx = [0.01 0.005]*100;
%     cp.u_MC_Rx = [0.01 0.005]*100;
%     cp.sigma_MC_Tx = 0.2;
%     cp.sigma_MC_Rx = 0.2;
    
    % PA
    cp.pa_gain = 25;

    % AC (array calibration)
    
    % AR (Array Response)
    cp.sigma_AR_Tx = 0.002;
    cp.sigma_AR_Rx = 0.002;
%     cp.sigma_AR_Tx = 0.2;
%     cp.sigma_AR_Rx = 0.2;
    cp.AR_B_Rx = eye(size(cp.H,1)); % could be cosine or omni-directional
    cp.AR_B_Tx = eye(size(cp.H,2));
    
    % AG (Array Geometry)
    cp.sigma_AG_Tx = 0.005;
    cp.sigma_AG_Rx = 0.005;
    % debug value
%     cp.sigma_AG_Tx = 0.05;
%     cp.sigma_AG_Rx = 0.05;

    % IQI
    cp.IQI_phi_Tx = 0.05;
    cp.IQI_m_Tx = 0.95;
    cp.IQI_phi_Rx = 0.05;
    cp.IQI_m_Rx = 0.95;
    
    
end