function c = default_THz_2D_SWM_parameters()

% signal    
    c.c = 3e8;
    c.fc = 140e9; % mmWave
    
%     c.BW = 1000e6;
    c.BW = 400e6;
%     c.BW = 200e6;

    c.K = 10;
    c.G = 5;

%     c.BW = 400e6;
%     c.K = 20;
    
    c.v0 = 0;   % speed
    c.AOD = 0;    % angle of arrival

    % Simulation parameters
    c.PU = [3 5]';
    
    c.D_Tx = 0;
    c.D_Rx = (0:31)';  % uniform linear array
    c.array_structure = "Digital";     % Digital
    c.wave_type = "PWM";
    
    % Model mismatch: SNS, SWM, BSE
    % spatial non-stationarity, spherical wave model, beam squint effect
    c.NF_SNS = "False";     % True, False
    c.NF_SWM = "False";
    c.NF_BSE = "False";
%     
%     c.NF_SNS = "True";     % True, False
%     c.NF_SWM = "True";
%     c.NF_BSE = "True";
    c.P = 10;
end