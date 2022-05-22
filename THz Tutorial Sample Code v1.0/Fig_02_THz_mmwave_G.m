% close all;
clear all;
clc;
% rng(0);
% sequence: k, b, m, r
%% Initialization
c = defaultParameters();
c.random_seed = "Fixed";
c.rng_ind = 1;
% BM & BRM channel
c.NB_dim = [4 4]';   % Dim of BS Elements
c.NR_dim = [4 4]';   
c.NM_dim = [4 4]';
c.PB = [0 0 0]';
c.PR = [5, 5, 0]';
c.PM = [10, 0, 0]';
c.OB = [0, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [150, 0, 0]';

c.PC = [];  % assign an empty matrix if not available
c.PR = [];  % assign an empty matrix if not available
c.G = 10;
% signal parameters
% mmWave
% c.fc = 60e9;
% c.BW = 100e6;    % 100M bandwidth
% c.K = 10;

% THz
c.fc = 300e9;
c.BW = 100e6;    % 100M bandwidth
% c.P = 1;   % 1 dBm
c.K = 10;

% SA: analog beamforming
c.NaeB = 1;     % antenna elements (AEs) per SA at BS
c.NaeM = 1;
c.BFmatB = zeros(2, prod(c.NB_dim))+[0 0]';  % row1: phi, row2: theta
c.BFmatM = zeros(2, prod(c.NM_dim))+[0 0]';
c.deltaDaeCoe = 1;    % in [lambda/2]. coe=1 --> 0.5*lambdac.


% Model parameters
c.symbol_type = "Random"; % Ones, Random, RandomSymbol, 
c.pos_type = "2D";         % 2D, 3D positions
c.ori_type = "1D";
c.syn_type = "Syn"; % Syn, Asyn

% c.AoSA_type = "SSP"; % SWM/PWM model for: SA phase/SA amplitude/AE
c.AoSA_type = "SPP"; % SWM/PWM model for: SA phase/SA amplitude/AE
c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized

% different realizations
c = updateParameters(c);


size(c.Xmk)
norm(c.Xmk, 'fro')
c0 = c;

%% SIM-1: SA size  (lambda/2, lambda/8)
% close all;
Iter_vec = 1:12;
PEB_cell = cell(1, length(Iter_vec));
OEB_cell = cell(1, length(Iter_vec));
sim_times = 100;
AEB_cell = cell(1, length(Iter_vec));
DEB_cell = cell(1, length(Iter_vec));

% c.symbol_type = "Random"; % Ones, Random, RandomSymbol
tic
parfor iter = 1:length(Iter_vec)
    disp("Iteration: " + num2str(iter));

    for sim_i = 1:sim_times
%         rng_ind = 0;
        if(iter == 16 || iter == 8)&& mod(sim_i, 10) == 0 
            disp("Simtimes: " + num2str(sim_i));
        end
        % Realization 1: mmWave benchmark
        
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.G = Iter_vec(iter);
        c.fc = 60e9;
        c.BW = 100e6;    % 100M bandwidth
        c.NB_dim = [4 4]';   % Dim of BS Elements
        c.NM_dim = [2 2]';
        sa_dim_B = 1;
        sa_dim_M = 1;
        c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [sa_dim_M sa_dim_M];
        c.DsaB = sa_dim_B;
        c.DsaM = sa_dim_M;
        c.BFmatB = zeros(2, prod(c.NB_dim));
        c.BFmatM = zeros(2, prod(c.NM_dim));
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(1, sim_i) = PEB;
        OEB_cell{iter}(1, sim_i) = OEB;
        
        % Realization 2: THz 2 2
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.G = Iter_vec(iter);
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        c.NB_dim = [4 4]';   % Dim of BS Elements
        c.NM_dim = [2 2]';
        sa_dim_B = 2;
        sa_dim_M = 2;
        c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [sa_dim_M sa_dim_M];
        c.DsaB = sa_dim_B;
        c.DsaM = sa_dim_M;
        c.BFmatB = zeros(2, prod(c.NB_dim));
        c.BFmatM = zeros(2, prod(c.NM_dim));
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(2, sim_i) = PEB;
        OEB_cell{iter}(2, sim_i) = OEB;

        % Realization 3: THz 5 5
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.G = Iter_vec(iter);
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        c.NB_dim = [4 4]';   % Dim of BS Elements
        c.NM_dim = [2 2]';
        sa_dim_B = 5;
        sa_dim_M = 5;
        c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [sa_dim_M sa_dim_M];
        c.DsaB = sa_dim_B;
        c.DsaM = sa_dim_M;
        c.BFmatB = zeros(2, prod(c.NB_dim));
        c.BFmatM = zeros(2, prod(c.NM_dim));
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(3, sim_i) = PEB;
        OEB_cell{iter}(3, sim_i) = OEB;
                
        % Realization 4: THz 10 10
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.G = Iter_vec(iter);
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        c.NB_dim = [4 4]';   % Dim of BS Elements
        c.NM_dim = [2 2]';
        sa_dim_B = 10;
        sa_dim_M = 10;
        c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [sa_dim_M sa_dim_M];
        c.DsaB = sa_dim_B;
        c.DsaM = sa_dim_M;
        c.BFmatB = zeros(2, prod(c.NB_dim));
        c.BFmatM = zeros(2, prod(c.NM_dim));
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(4, sim_i) = PEB;
        OEB_cell{iter}(4, sim_i) = OEB;

        % Realization 5: THz, Asyn, AOSA
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.G = Iter_vec(iter);
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        c.K = 10;
        c.NB_dim = [4 4]';   % Dim of BS Elements
        c.NM_dim = [1 1]';
        sa_dim_B = 5;
        sa_dim_M = 5;
        c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [sa_dim_M sa_dim_M];
        c.DsaB = sa_dim_B;
        c.DsaM = sa_dim_M;
        c.BFmatB = zeros(2, prod(c.NB_dim));
        c.BFmatM = zeros(2, prod(c.NM_dim));
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(5, sim_i) = PEB;
        OEB_cell{iter}(5, sim_i) = OEB;
    end
end
toc

% save data_Fig_2.mat

% load data_Fig_2.mat

%% plot
PEB_data = zeros(5, length(Iter_vec));
OEB_data = zeros(5, length(Iter_vec));

for i = 1:length(Iter_vec)
    PEB_data(:, i) = mean(abs(PEB_cell{i}),2);
    OEB_data(:, i) = mean(abs(OEB_cell{i}),2);
end


linewidth = 1.5;
markersize = 8;
% % PEB
figure;semilogy(Iter_vec, PEB_data(1,:), 'k-', 'Linewidth', 1.5);
hold on;semilogy(Iter_vec, PEB_data(2,:), 'rs--', 'Linewidth', linewidth, 'MarkerSize', markersize);
hold on;semilogy(Iter_vec, PEB_data(5,:), 'b^--', 'Linewidth', linewidth, 'MarkerSize', markersize);
hold on;semilogy(Iter_vec, PEB_data(3,:), 'bo-', 'Linewidth', linewidth, 'MarkerSize', markersize);
hold on;semilogy(Iter_vec, PEB_data(4,:), 'gd--', 'Linewidth', linewidth, 'MarkerSize', markersize);
% hold on;semilogy(Iter_vec, PEB_data(6,:), 'b--', 'Linewidth', linewidth);

legend('{mmWave ($N_B=4\times4$, $N_U=2\times2$)}', ...
    '{THz ($N_B=4\times4$, $N_U=2\times2$, $\mathring N_B=\mathring N_U$=2$\times$2)}', ...
    '{THz ($N_B=4\times4$, $N_U=1\times1$, $\mathring N_B=\mathring N_U$=5$\times$5)}', ...
    '{THz ($N_B=4\times4$, $N_U=2\times2$, $\mathring N_B=\mathring N_U$=5$\times$5)}', ...
    '{THz ($N_B=4\times4$, $N_U=2\times2$, $\mathring N_B=\mathring N_U$=10$\times$10)}', ...
    'Location', 'NorthEast', 'Interpreter','Latex');
   
xlabel('{Number of Transmissions $\mathcal{G}$}','Interpreter','Latex')
ylabel('PEB [m]','Interpreter','Latex');
axis([1 max(Iter_vec) 0.07 4]);
grid on;
set(gca,'fontsize', 16);
set(gcf,'position', [100,100,410*1.4, 410])

% print -dpng -r600 sim_2_transmissions_PEB.png

% saveas(gcf,'sim_2_transmissions_PEB', 'epsc')

%% OEB [^\circ]
figure;semilogy(Iter_vec, OEB_data(1,:), 'k-', 'Linewidth', 1.5);
hold on;semilogy(Iter_vec, OEB_data(5,:), 'rs-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(2,:), 'rs--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(3,:), 'bo-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(4,:), 'gd--', 'Linewidth', linewidth);


% legend( '{60\hspace{0.5mm}GHz/0.1GHz/$\bar{N}_A$=1}', ...
%         '{0.3THz/0.1GHz/$\bar{N}_A$=1}', '{0.3THz/0.1GHz/$\bar{N}_A$=1\hspace{1.7mm}/$D_{SA}$=2.5$\lambda$}', ...
%         '{0.3THz/0.1GHz/$\bar{N}_A$=25/$D_{SA}$=2.5$\lambda$}', ...
%         '{0.3THz/\hspace{2.7mm}1GHz/$\bar{N}_A$=$25$/$D_{SA}$=2.5$\lambda$}',...
%         'Location', 'NorthEast', 'Interpreter','Latex');
   
% xlabel('\textbf{Number of BS Antennas}','Interpreter','Latex')
% ylabel('\textbf{OEB [deg]}','Interpreter','Latex');
xlabel('{Number of Transmissions $\mathcal{G}$}','Interpreter','Latex')
ylabel('{OEB [deg]}','Interpreter','Latex');
% axis([0 100 10^-2 10^2]);

grid on;
set(gca,'fontsize', 14);
set(gcf,'position', [100,100, 400*1.37, 400])

% print -dpng -r600 sim_2_transmissions_OEB.png

% saveas(gcf,'sim_2_transmissions_OEB', 'epsc')

%%


