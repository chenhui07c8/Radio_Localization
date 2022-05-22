% Matlab Code for Paper Published in IEEE Communication Surveys & Tutorials
% "A Tutorial on Terahertz-Band Localization for 6G Communication Systems"
% Author: Hui Chen
% Email: hui.chen@chalmers.se; hui.chen@kaust.edu.sa

close all;
clear all;
clc;
folder = 'THz Loc CRLB_v1.0';
addpath(genpath(folder));


%% Initialization
c = defaultParameters();
c.rng_ind = 0;
c.random_seed = "Fixed";   % "Random", "Fixed"

% BM & BRM channel
c.NB_dim = [4 4]';   % Dim of BS Elements
c.NR_dim = [4 4]';   
c.NM_dim = [2 2]';
c.PB = [0 0 0]';
c.PR = [5, 5, 0]';
c.PM = [10, 0, 0]';
c.OB = [0, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [150, 0, 0]';

c.PC = [];  % assign an empty matrix if not available
c.PR = [];  % assign an empty matrix if not available
c.G = 10;   % transmission times

% signal parameters
% mmWave
c.fc = 60e9;
c.BW = 100e6;    % 100M bandwidth
c.K = 10;
c.NsaB_dim = [1 1];     % antenna elements (AEs) per SA at BS
c.NsaM_dim = [1 1];

% THz
% c.fc = 300e9;
% c.BW = 100e6;    % 100M bandwidth
% c.P = 1;   % 1 dBm

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
%% RIS coefficients & channel realization
% clc
uB = c.uB;
uBM = c.uBM;
uBRM = c.uBRM;
uBCM = c.uBCM;

disp(['Value BM : ' num2str(norm(uBM, 'fro'))])
disp(['Value BRM: ' num2str(norm(uBRM, 'fro'))])
disp(['Value BCM: ' num2str(norm(uBCM, 'fro'))])
disp(['Value SUM: ' num2str(norm(uBM + uBRM + uBCM, 'fro'))])

efficiency = norm(uBM + uBRM, 'fro')/(norm(uBM, 'fro') + norm(uBRM, 'fro'));
disp(['Improvement: ' num2str(efficiency)]);

disp('Value SUM should have the largest value, hopefully.')

%%
c = c0;
c.fc = 60e9;
c.BW = 100e6;    % 100M bandwidth
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[~, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));

%% Simulation Start
% close all;
Iter_vec = (1:10).^2;
PEB_cell = cell(1, length(Iter_vec));
OEB_cell = cell(1, length(Iter_vec));
sim_times = 10;
AEB_cell = cell(1, length(Iter_vec));
DEB_cell = cell(1, length(Iter_vec));

tic
parfor iter = 1:length(Iter_vec)

    for sim_i = 1:sim_times
        disp("Iteration: " + num2str(iter) + " | Simulations:" + num2str(sim_i));
        
        % Realization 1: mmWave, Asyn
        c = c0;
        c.NB_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';
        c.fc = 60e9;
        c.BW = 100e6;    % 100M bandwidth
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [~, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(1, sim_i) = PEB;
        OEB_cell{iter}(1, sim_i) = OEB;
        CRLB2 = FIM(1:end-1, 1:end-1)^(-1);
        AEB = sqrt(CRLB2(3,3))/pi*180;
        DEB = sqrt(CRLB2(5,5))*3e8;
        AEB_cell{iter}(1, sim_i) = AEB;
        DEB_cell{iter}(1, sim_i) = DEB;
        
        % Realization 2: THz, Asyn
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.NB_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [~, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(2, sim_i) = PEB;
        OEB_cell{iter}(2, sim_i) = OEB;
        CRLB2 = FIM(1:end-1, 1:end-1)^(-1);
        AEB = sqrt(CRLB2(3,3))/pi*180;
        DEB = sqrt(CRLB2(5,5))*3e8;
        AEB_cell{iter}(2, sim_i) = AEB;
        DEB_cell{iter}(2, sim_i) = DEB;
        
        % Realization 3: THz, Asyn, AOSA
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.NB_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        sa_dim_B = 1;
        sa_dim_M = 1;
        c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [sa_dim_M sa_dim_M];
        c.DsaB = 5;
        c.DsaM = 5;
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [~, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(3, sim_i) = PEB;
        OEB_cell{iter}(3, sim_i) = OEB;
        CRLB2 = FIM(1:end-1, 1:end-1)^(-1);
        AEB = sqrt(CRLB2(3,3))/pi*180;
        DEB = sqrt(CRLB2(5,5))*3e8;
        AEB_cell{iter}(3, sim_i) = AEB;
        DEB_cell{iter}(3, sim_i) = DEB;
 
        % Realization 4: THz, Syn, AOSA, BF
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.NB_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';
        c.fc = 300e9;
        c.BW = 100e6;    % 100M bandwidth
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];
        c.DsaB = 5;
        c.DsaM = 5;
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [~, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(4, sim_i) = PEB;
        OEB_cell{iter}(4, sim_i) = OEB;
        CRLB2 = FIM(1:end-1, 1:end-1)^(-1);
        AEB = sqrt(CRLB2(3,3))/pi*180;
        DEB = sqrt(CRLB2(5,5))*3e8;
        AEB_cell{iter}(4, sim_i) = AEB;
        DEB_cell{iter}(4, sim_i) = DEB;
        
        % Realization 5: THz, Syn, AOSA, BW
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.NB_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';
        c.K = 100;
        c.fc = 300e9;
        c.BW = 1000e6;    % 100M bandwidth
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [~, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(5, sim_i) = PEB;
        OEB_cell{iter}(5, sim_i) = OEB;
        CRLB2 = FIM(1:end-1, 1:end-1)^(-1);
        AEB = sqrt(CRLB2(3,3))/pi*180;
        DEB = sqrt(CRLB2(5,5))*3e8;
        AEB_cell{iter}(5, sim_i) = AEB;
        DEB_cell{iter}(5, sim_i) = DEB;

        % Realization 6: Beamforming to MD
        c = c0;
        c.rng_ind = iter*sim_times + sim_i;
        c.NB_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';
        c.K = 100;
        c.fc = 300e9;
        c.BW = 1000e6;    % 100M bandwidth
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c = updateParameters(c);
        cM = get_multiple_transmission_parameters(c);
        cM.BFmatB(1,:,:) = c.phiBM_loc;
        cM.BFmatM(1,:,:) = c.phiMB_loc;
        cM.BFmatB(1,:,:) = (rand(c.NB, c.G)-0.5)*2*10;
        cM.BFmatM(1,:,:) = c.phiMB_loc + (rand(c.NM, c.G)-0.5)*2*10;
        [FIM, ~] = get_FIM_multi_transmission(c, cM);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(6, sim_i) = PEB;
        OEB_cell{iter}(6, sim_i) = OEB;
        CRLB2 = FIM(1:end-1, 1:end-1)^(-1);
        AEB = sqrt(CRLB2(3,3))/pi*180;
        DEB = sqrt(CRLB2(5,5))*3e8;
        AEB_cell{iter}(6, sim_i) = AEB;
        DEB_cell{iter}(6, sim_i) = DEB;
    end
end
toc

% save data_Fig_1.mat

% load data_Fig_1.mat

%% plot
PEB_data = zeros(6, length(Iter_vec));
OEB_data = zeros(6, length(Iter_vec));

for i = 1:length(Iter_vec)
    PEB_data(:, i) = mean(abs(PEB_cell{i}),2);
    OEB_data(:, i) = mean(abs(OEB_cell{i}),2);
end

linewidth = 1.5;
% PEB
figure;semilogy(Iter_vec, PEB_data(1,:), 'k-', 'Linewidth', 1.5);
hold on;semilogy(Iter_vec, PEB_data(2,:), 'r+--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, PEB_data(3,:), 'ro--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, PEB_data(4,:), 'b+-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, PEB_data(5,:), 'bo-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, PEB_data(6,:), 'b--', 'Linewidth', linewidth);

% 1/mean(PEB_data(5,2:end)./PEB_data(1,2:end))
% 1/mean(PEB_data(6,2:end)./PEB_data(1,2:end))
% 4.88, 26.95

% 1/mean(OEB_data(5,2:end)./OEB_data(1,2:end))
% 1/mean(OEB_data(6,2:end)./OEB_data(1,2:end))
% 6.01
% 24.7

legend( '{mmW MIMO ($f_c$=60GHz)}', ...
    '{THz-1 (mmW setup with $f_c$=0.3THz)}', ...
    '{THz-2 (THz-1 with $\Delta$=2.5$\lambda$)}', ...
    '{THz-3 (THz-2 with AOSA, $\mathring{N}_A$=5$\times5$)}', ...
    '{THz-4 (THz-3 with $W$=1GHz)}',...
    '{THz-5 (THz-4 with Prior Info)}',...
    'Location', 'NorthEast', 'Interpreter','Latex');
   
xlabel('{Number of Antennas/SAs at BS ($N_\mathrm{B}$)}','Interpreter','Latex')
ylabel('PEB [m]','Interpreter','Latex');
axis([0 100 10^-2.6 10^2]);
grid on;
set(gca,'fontsize', 18);
set(gcf,'position', [100,100,400*1.4, 400])


% OEB [^\circ]
figure;semilogy(Iter_vec, OEB_data(1,:), 'k-', 'Linewidth', 1.5);
hold on;semilogy(Iter_vec, OEB_data(2,:), 'r+--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(3,:), 'ro--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(4,:), 'b+-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(5,:), 'bo-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, OEB_data(6,:), 'b--', 'Linewidth', linewidth);

% legend( '{60\hspace{0.5mm}GHz/0.1GHz/$\bar{N}_A$=1}', ...
%         '{0.3THz/0.1GHz/$\bar{N}_A$=1}', '{0.3THz/0.1GHz/$\bar{N}_A$=1\hspace{1.7mm}/$D_{SA}$=2.5$\lambda$}', ...
%         '{0.3THz/0.1GHz/$\bar{N}_A$=25/$D_{SA}$=2.5$\lambda$}', ...
%         '{0.3THz/\hspace{2.7mm}1GHz/$\bar{N}_A$=$25$/$D_{SA}$=2.5$\lambda$}',...
%         'Location', 'NorthEast', 'Interpreter','Latex');
   
xlabel('{Number of Antennas/SAs at BS ($N_\mathrm{B})$}','Interpreter','Latex')
ylabel('{OEB [deg]}','Interpreter','Latex');
axis([0 100 10^-1.2 10^2]);

grid on;
set(gca,'fontsize', 18);
set(gcf,'position', [100,100, 400*1.4, 400])

% print -dpng -r600 sim_1_THz_vs_mmW_PEB.png

% print -dpng -r600 sim_1_THz_vs_mmW_OEB.png

% saveas(gcf,'sim_1_THz_vs_mmW_PEB', 'epsc')

% saveas(gcf,'sim_1_THz_vs_mmW_OEB', 'epsc')



