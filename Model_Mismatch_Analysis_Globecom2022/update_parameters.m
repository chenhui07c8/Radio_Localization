function c = update_parameters(c)

    c.delta_f = c.BW/c.K;
    % c.Tsym = 0.01;
    c.lambdac = c.c/c.fc;
    c.fk = c.fc + (0:c.K-1)*c.delta_f;
    c.lambdak = c.c./c.fk;
    c.L = 1;
    c.N = length(c.D_Rx);
    
        
    c.rho = c.lambdac/4/pi/norm(c.PU);
%     c.xi = rand(1,1)*2*pi;
    c.xi = 0;
    c.gain = c.rho*exp(1j*c.xi);
    
%     cp = get_Rx_symbol_PWM(cp);
%     cp = get_Rx_symbol_SWM(cp);

% geometry
    c.AOA = atan2d(c.PU(2), c.PU(1));    % angle of departure
    c.d = norm(c.PU);
    c.AOD = 0;

    
    c.D_Rx = c.D_Rx - mean(c.D_Rx);
    c.D_Rx_Loc = [zeros(1, length(c.D_Rx)); (c.D_Rx-mean(c.D_Rx))'*c.lambdac/2];
    
    c.N_Tx = length(c.D_Tx);
    c.N_Rx = length(c.D_Rx);
    
    rng(0);
    if c.array_structure == "Digital"
        c.M_Tx = c.N_Tx;
        c.M_Rx = c.N_Rx;
        c.WRx = ones(1, c.G);
        c.WTx = ones(1, c.G);
    elseif c.array_structure == "Analog"
        
        c.M_Tx = 1;
        c.M_Rx = 1;
        % frequency-independent phase-shifters
        c.WRx = exp(1j*rand(c.N_Rx, c.G)*2*pi)/sqrt(c.N_Rx); % same for all the subcarriers
        % frequency-dependent phase-shifters
        c.WRx = exp(1j*rand(c.N_Rx, c.G)*2*pi)/sqrt(c.N_Rx); % same for all the subcarriers
        
%         c.WTx = exp(1j*rand(c.N_Tx, c.G)*2*pi)/sqrt(c.N_Tx); % same for all the subcarriers
        c.WTx = ones(1, c.G);
    end
    
    c.s = ones(c.K, 1);    % transmitted symbols
%     c.s = exp(1j*2*pi*rand(c.K, 1))/sqrt(c.K);    % transmitted symbols
%     c.s = fft(c.s)/sqrt(c.K);    % transmitted symbols
%     norm(c.s)
%     c.X0 = c.s;

% hardware

    
% environment
%     c.fc = 100e9; % sub-THz
%     c.BW = 1000e6;
%     c.K = 100;
    
    % c.X0 = exp(1j*2*pi*rand(c.N_Tx, c.G))/sqrt(c.K);

   

    c.operationBW = c.BW;    % Operation bandwidth for Thermal noise
    c.K_boltzmann = 1.3806e-23;      % boltzmann constant
    c.Temperature = 298.15;         % temprature 25 celsius
    c.Pn = c.K_boltzmann*c.Temperature*c.operationBW*1000;    % thermal noise linear (in mW)
    c.Pn_dBm = 10*log10(c.Pn);      % thermal noise decibel (in dB)
    c.sigma_in = sqrt(c.Pn);          % Johnson–Nyquist noise: sigma^2 = N_0
    c.NoiseFigure = 10;       % noise figure 3dB.
    c.sigma = sqrt(10^(c.NoiseFigure/10))*c.sigma_in;

    
    
    
    
    c.snr = pow2db(1/c.sigma^2);  % SNR for P = 1
%     c.pa_gain = 25;

end