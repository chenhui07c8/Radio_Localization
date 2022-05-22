% close all;
clear all;
clc;
% rng(0);
% sequence naming priority: k, b, m, r
%% Initialization
c = defaultParameters();

c.rng_ind = 0;

c.NB_dim = [4 4]';   % Dim of BS Elements
c.NR_dim = [1 1]';   
c.NM_dim = [4 4]';
c.PB = [0 0 0]';
c.PR = 0.5*[5, 5, 0]';
c.PM = [2, 0, 0]';
c.OB = [60, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [180, 0, 0]';
c.P = 10;
c.K = 10;
c.PC = [];
c.PR = [];
% mmwave
% c.fc = 60e9;    % carrier frequency
% c.BW = 1000e6;    % 100M bandwidth

% THz
c.fc = 300e9;    % carrier frequency
c.BW = 10e9;    % 100M bandwidth
% c.K = 4;

sa_dim_B = 5;
sa_dim_M = 5;
% sa_dim_B = 4;
% sa_dim_M = 4;
c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
c.NsaM_dim = [sa_dim_M sa_dim_M];
c.DsaB = 5;
c.DsaM = 5;

% Hardware Parameters
% Antenna gain
% c.GB = [1 pi/2 pi/2];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;
% c.GR = [1 pi/3 pi/3];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;
% c.GM = [1 pi/2 pi/2];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;

% c.BFmatM = zeros(2, prod(c.NM_dim))+[0 0]';
% reshape([1 2 3; 4 5 6; 7 8 9]', [1, 9])

c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
[c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
[c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
% c.BFmatB(:,1) = [-5 0]';  % row1: phi, row2: theta
c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
[c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
% c.BFmatB(:,1) = [-5 0]';  % row1: phi, row2: theta


% Model Parameters
c.symbol_type = "Random"; % Ones, Random, RandomSymbol
% c.BeamSplit = "True";
% c.signal_model = "PWM";
c.pos_type = "2D";         % 2D, 3D positions
c.ori_type = "1D";
c.syn_type = "Asyn"; % Syn, Asyn
% c.syn_type = "Syn"; % Syn, Asyn

% c.AoSA_type = "SSP"; % SWM/PWM model for: SA phase/SA amplitude/AE
c.AoSA_type = "SPP"; % SWM/PWM model for: SA phase/SA amplitude/AE

c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
c.fixed_noise = "False";         % False, True

    
c = updateParameters(c);
c.lambdac

size(c.Xmk)
% c.Xmn
norm(c.Xmk, 'fro')

c0 = c;
% c.GantBM
%% RIS
uB = c.uB;
uBM = c.uBM;
uBRM = c.uBRM;
uBCM = c.uBCM;
%
n = 1
b = 1;
disp(['Phase BM : ' num2str(phase(uBM(n,b)))]);
disp(['Phase BRM: ' num2str(phase(uBRM(n,b)))]);

disp(['Value BM : ' num2str(norm(uBM(n, :)))])
disp(['Value BRM: ' num2str(norm(uBRM(n, :)))])
disp(['Value SUM: ' num2str(norm(uBM(n, :) + uBRM(n, :)))])


disp(['Value BM : ' num2str(norm(uBM, 'fro'))])
disp(['Value BRM: ' num2str(norm(uBRM, 'fro'))])
disp(['Value BCM: ' num2str(norm(uBCM, 'fro'))])
disp(['Value SUM: ' num2str(norm(uBM + uBRM, 'fro'))])


efficiency = norm(uBM + uBRM, 'fro')/(norm(uBM, 'fro') + norm(uBRM, 'fro'));
disp(['Improvement: ' num2str(efficiency)]);

disp('Value SUM should have the largest value, hopefully.')

% energyB0 = norm(uBM, 'fro');
% energyR0 = norm(uBRM, 'fro');
uB(end,end)

%% CRLB
% c.pos_type = "2D";          % 2D, 3D positions
% c.ori_type = "1D";          % 1D, 2D, 3D orientations
c = c0;
c.beam_split = "False";          % True, False
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));


c = c0;
c.beam_split = "True";          % True, False
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));

min(c.AeqBM(:))
max(c.AeqBM(:))


%% CRLB
Iter_vec = logspace(-1, log10(10), 10);
Iter_vec = round(Iter_vec*100)/100;
sim_times = 100;

PEB_cell = cell(1, length(Iter_vec));
OEB_cell = cell(1, length(Iter_vec));


parfor iter = 1:length(Iter_vec)
    disp("Iteration: " + num2str(iter));
    PEB_cell{iter}(1:6, 1:sim_times) = 0;
    OEB_cell{iter}(1:6, 1:sim_times) = 0;
    for sim_i = 1:sim_times

        % Realization 1: benchmark without NLOS
        c = c0;
        c.rng_ind = sim_i;
        c.PC = [];
        c.BW = Iter_vec(iter)*1e9;
        c.OB = [45, 0, 0]';
        c.K = c.BW/(10e6);
        c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
        c.beam_split = "True";          % True, False
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 5;
        c.DsaM = 5;
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(1,sim_i) = PEB;
        OEB_cell{iter}(1,sim_i) = OEB;

        % Realization 2: benchmark without NLOS
        c = c0;
        c.rng_ind = sim_i;
        c.PC = [];
        c.BW = Iter_vec(iter)*1e9;
        c.OB = [45, 0, 0]';
        c.K = c.BW/(10e6);
        c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 5;
        c.DsaM = 5;
        c.beam_split = "False";          % True, False
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(2,sim_i) = PEB;
        OEB_cell{iter}(2,sim_i) = OEB;
        
        
        % Realization 2: benchmark without NLOS
        c = c0;
        c.rng_ind = sim_i;
        c.PC = [];
        c.BW = Iter_vec(iter)*1e9;
        c.OB = [15, 0, 0]';
        c.K = c.BW/(10e6);
        c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 5;
        c.DsaM = 5;
        c.beam_split = "True";          % True, False
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(3,sim_i) = PEB;
        OEB_cell{iter}(3,sim_i) = OEB;
        
        
        
        % Realization 2: benchmark without NLOS
        c = c0;
        c.rng_ind = sim_i;
        c.PC = [];
        c.BW = Iter_vec(iter)*1e9;
        c.OB = [15, 0, 0]';
        c.K = c.BW/(10e6);
        c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
        
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 5;
        c.DsaM = 5;
        c.beam_split = "False";          % True, False
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(4,sim_i) = PEB;
        OEB_cell{iter}(4,sim_i) = OEB;
        
        
        % Realization 2: benchmark without NLOS
        c = c0;
        c.rng_ind = sim_i;
        c.PC = [];
        c.BW = Iter_vec(iter)*1e9;
        c.OB = [15, 0, 0]';
        c.K = c.BW/(10e6);
        c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
        
        c.NsaB_dim = [2 2];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 2;
        c.DsaM = 5;        
        c.beam_split = "True";          % True, False
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(5,sim_i) = PEB;
        OEB_cell{iter}(5,sim_i) = OEB;
        
               
        % Realization 2: benchmark without NLOS
        c = c0;
        c.rng_ind = sim_i;
        c.PC = [];
        c.BW = Iter_vec(iter)*1e9;
        c.OB = [15, 0, 0]';
        c.K = c.BW/(10e6);

        c.tBM = (c.PM-c.PB)/norm(c.PM-c.PB);
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
        
        c.NsaB_dim = [2 2];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 2;
        c.DsaM = 5;        
        c.beam_split = "False";          % True, False
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(6,sim_i) = PEB;
        OEB_cell{iter}(6,sim_i) = OEB;
    end
end

% save data_Fig_4.mat

% load data_Fig_4.mat

%% plot
PEB_data = zeros(6, length(Iter_vec));
OEB_data = zeros(6, length(Iter_vec));

for i = 1:length(Iter_vec)
    PEB_data(:, i) = mean(PEB_cell{i},2);
    OEB_data(:, i) = mean(OEB_cell{i},2);
end

% for i = 1:length(Iter_vec)
%     PEB_data(:, i) = PEB_cell{i};
%     OEB_data(:, i) = OEB_cell{i};
% end

linewidth = 1.5;
% PEB
figure;loglog(Iter_vec, PEB_data(5,:), 'g+-', 'Linewidth', linewidth);
hold on;loglog(Iter_vec, PEB_data(6,:), 'g^--', 'Linewidth', linewidth);
hold on;loglog(Iter_vec, PEB_data(1,:), 'b+-', 'Linewidth', linewidth);
hold on;loglog(Iter_vec, PEB_data(2,:), 'bo--', 'Linewidth', linewidth);
hold on;loglog(Iter_vec, PEB_data(3,:), 'r+-', 'Linewidth', linewidth);
hold on;loglog(Iter_vec, PEB_data(4,:), 'rs--', 'Linewidth', linewidth);

legend('{$\mathring{N}_\mathrm{B}$= 4\hspace{0.8mm}/$\alpha_B$=15$^\circ$/with BSE}', ...
       '{$\mathring{N}_\mathrm{B}$= 4\hspace{0.8mm}/$\alpha_B$=15$^\circ$/w/o BSE}', ...
       '{$\mathring{N}_\mathrm{B}$=25/$\alpha_B$=45$^\circ$/with BSE}', ...
       '{$\mathring{N}_\mathrm{B}$=25/$\alpha_B$=45$^\circ$/w/o BSE}', ...
       '{$\mathring{N}_\mathrm{B}$=25/$\alpha_B$=15$^\circ$/with BSE}',...
       '{$\mathring{N}_\mathrm{B}$=25/$\alpha_B$=15$^\circ$/w/o BSE}', ...
       'Location', 'SouthWest', 'Interpreter','Latex');


xlabel('{Bandwidth [GHz]}','Interpreter','Latex')
ylabel('PEB [m]','Interpreter','Latex');

axis([min(Iter_vec) max(Iter_vec) 10^-2.4 10^0]);
grid on;
set(gca,'fontsize', 16);
set(gcf,'position', [100,100,350*1.43, 350])

% print -dpng -r600 sim_4_beamsplit_PEB.png

% saveas(gcf,'sim_4_beamsplit_PEB', 'epsc')


