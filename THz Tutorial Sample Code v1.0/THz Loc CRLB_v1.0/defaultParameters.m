
% Initialize default channel parameters, including several parts.
% 1). Model Parameters: set channel and optimization features
% 2). Geometry Parameters: set position and orientation
% 3). Infrastructure Parameters: array & AOSA dimensions, signal and antenna parameters
% 4). Environment Parameters: noise figure, noise level


% channel parameters as a structure 'c'.
function c = defaultParameters()

% Model Parameters
    c.pos_type = "2D";          % 2D, 3D positions
    c.ori_type = "1D";          % 1D, 2D, 3D orientations
    c.channel_knowledge = "False";  % True , False, Partial
    c.syn_type = "Syn";             % Syn, Asyn
    c.fixed_noise = "False";        % True for user defined. False for operation BW.
    c.symbol_type = "Random";         % Ones, Random, RandomSymbol
    c.beam_split = "True";          % True, False
    c.wave_model = "SWM";       % PWM, SWM
    c.multipleTransmission = 'RRR';  % Beamforming angles for multiple transmission
    c.AoSA_type = "SSP"; % SWM/PWM model for: SA phase/SA amplitude/AE
    c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Fixed, Customized
    c.RIS_quantization = 0;
    c.gain_type = "Localization";  % 'Localization': one gain one array channel; 'Accurate': one gain one SA channel.
    c.random_seed = "Random";   % "Random", "Fixed"
    
% Geometry Parameters
    % global position: 3x1 vector (2D by setting P(3) as zero)
    c.PB = [0 0 0]';
    c.PM = [4, 0, 0]';
    c.PC = [];
    c.PR = [2 2]';
    % rotation matrix: 3x1 vector (1D by setting P(2) & P(3) as zeros)
    c.OB = [0, 0, 0]';
    c.OR = [-90, 0, 0]';
    c.OM = [180, 0, 0]';
    
% Infrastructure Parameters
    c.lightspeed = 3e8;      % speed of light
    % Array Dimensions
    c.NB_dim = [3 3]';   % [# of rows, # of columns]
    c.NR_dim = [3 3]';
    c.NM_dim = [3 3]';

    c.NB = c.NB_dim(1)*c.NB_dim(2);       % # of BS Elements
    c.NR = c.NR_dim(1)*c.NR_dim(2);       % # of RIS Elements
    c.NM = c.NM_dim(1)*c.NM_dim(2);       % # of MD Elements

    c.DsaR = 1;     % interval of RIS elements (in lambda/2).
    c.DsaB = 1;
    c.DsaM = 1;

    % AOSA: assuming a planar array
    c.NsaB_dim = [1 1];     % Dimension of SA at BS
    c.NsaM_dim = [1 1]; 
    c.NsaR_dim = [1 1];
    
    c.NsaB = prod(c.NsaB_dim);  % number of AE per SA at BS
    c.NsaM = prod(c.NsaM_dim);
    c.NsaR = prod(c.NsaR_dim);
    c.DaeB = 1;	% Element spacing coefficients (BS/RIS/M)
    c.DaeR = 1;	% AE spacing coefficients
    c.DaeM = 1;	% AE spacing coefficients
    
    c.BFmatB = zeros(2, c.NB)+[0 0]';  % row1: phi, row2: theta
    c.BFmatM = zeros(2, c.NM)+[0 0]';
    c.BFmatRB = zeros(2, c.NR)+[0 0]';  % analoy beamforming angle of the input signal
    c.BFmatRM = zeros(2, c.NR)+[0 0]';  % analog beamforming angle of the output signal
    
    % Signal & Channel Parameters
    c.fc = 60e9;    % carrier frequency
    c.lambdac = c.lightspeed/c.fc;
    c.BW = 100e6;   % 100M bandwidth
    c.Kabs = 16;    % atmosphere attenuation 16 dB/km
    c.P = 10;        % in [mW]. 1000mW = 30dBm. 1mW = 0dBm
    c.beta = 2;     % symchronization offset (in [meter])
    c.K = 4;        % # of carrier frequencies

    % Antenna gain. [gain, max_phi, max_theta], HPBW_phi = 2*max_phi.
    c.GB = [1 180 180];  % [1 pi pi] for omnidirectional antennas
    c.GM = [1 180 180];
    c.GR = [1 180 180];
    
%     c.GB = [1 30 30]; % direcitonal antennas with HPBW = 60
%     c.GM = [1 30 30];

% Environment Parameters
    c.operationBW = 1e8;    % Operation bandwidth for Thermal noise
    c.Kb = 1.3806e-23;      % boltzmann constant
    c.Temperature = 298.15;         % temprature 25 celsius
    c.Pn = c.Kb*c.Temperature*c.operationBW*1000;    % thermal noise linear (in mW)
    c.Pn_dBm = 10*log10(c.Pn);      % thermal noise decibel (in dB)
    c.sigma_in = sqrt(c.Pn);          % Johnson–Nyquist noise: sigma^2 = N_0
    c.NoiseFigure = 13;       % noise figure 3dB.
    c.sigma = sqrt(10^(c.NoiseFigure/10))*c.sigma_in;      % Output noise level
    % c.P/c.sigma 


end
