% Matlab Code for Paper Published in IEEE Communication Surveys & Tutorials
% "A Tutorial on Terahertz-Band Localization for 6G Communication Systems"
% Author: Hui Chen
% Email: hui.chen@chalmers.se; hui.chen@kaust.edu.sa

% close all;
clear all;
clc;
folder = 'THz Loc CRLB_v1.0';
addpath(genpath(folder));

% rng(0);
% sequence: k, b, m, r

%% Initialization
c = defaultParameters();

% BM & BRM channel
c.NB_dim = [4 4]';   % Dim of BS Elements
c.NR_dim = [4 4]';   
c.NM_dim = [4 4]';
c.PB = [0 0 0]';
c.PR = 0.1*[5, 5, 1]';
c.PM = 0.1*[5, 4, 0.5]';


c.OB = [0, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [150, 0, 0]';
c.P = 1;
c.rng_ind = 1;
c.G = 10;
% BCM channels; NLOS/clusters
% c.PC = [4 -5 0; 6 -5 0; 2 7 0]';
% c.coeClu = 10*[0.2 0.2 0.5];
% c.PC = [5 -3.000 0]';
% c.coeClu = 2*[0.2];
c.PC = [];  % assign an empty matrix if not available

% signal parameters
% mmWave
c.fc = 60e9;
c.BW = 100e6;    % 100M bandwidth
c.K = 10;

% THz
% c.fc = 300e9;
% c.BW = 100e6;    % 100M bandwidth
% c.K = 10;



% SA: analog beamforming
sa_dim_B = 1;     % antenna elements (AEs) per SA at BS
sa_dim_M = 1;
c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
c.NsaM_dim = [sa_dim_M sa_dim_M];
c.BFmatB = zeros(2, prod(c.NB_dim))+[0 0]';  % row1: phi, row2: theta
c.BFmatM = zeros(2, prod(c.NM_dim))+[0 0]';


% Model parameters
c.symbol_type = "Random"; % Ones, Random, RandomSymbol
% c.BeamSplit = "True";
% c.signal_model = "PWM";
c.pos_type = "2D";         % 2D, 3D positions
c.ori_type = "1D";
% c.pos_type = "2D";         % 2D, 3D positions
% c.ori_type = "1D";
c.syn_type = "Asyn"; % Syn, Asyn

c.AoSA_type = "SSP"; % SWM/PWM model for: SA phase/SA amplitude/AE
% c.AoSA_type = "SPP"; % SWM/PWM model for: SA phase/SA amplitude/AE

c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
% c.RIS_optimizer = "Random";   % Elzanaty, He, Random, Quantized, Fixed, Customized
c.fixed_noise = "False";         % false, true
% c.operationBW = 10e6;

% different realizations
c = updateParameters(c);
% c.BFmatM(1,:) = c.phiMB_loc;


size(c.Xmk)
norm(c.Xmk, 'fro')
c0 = c;
% RIS coefficients & channel realization
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

%% CRLB
% c.pos_type = "2D";          % 2D, 3D positions
% c.ori_type = "1D";          % 1D, 2D, 3D orientations

c = c0;
c.fc = 60e9;
c.BW = 100e6;    % 100M bandwidth
c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));

c = c0;
c.fc = 60e9;
c.BW = 100e6;    % 100M bandwidth
c.RIS_optimizer = "He";   % Elzanaty, He, Random, Quantized, Fixed, Customized
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));

%% SIM-1: SA size (lambda/2, lambda/8)
% close all;
sim_times = 1;
Iter_vec = (1:20).^2;

PEB_cell = cell(1, length(Iter_vec));
OEB_cell = cell(1, length(Iter_vec));

parfor iter = 1:length(Iter_vec)
    disp("Iteration: " + num2str(iter));
    for sim_i = 1:sim_times

        % Realization 1: mmWave
        c = c0;
        c.rng_ind = sim_i;
        c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
        c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(1, sim_i) = PEB;
        OEB_cell{iter}(1, sim_i) = OEB;
        
        % Realization 2: 16 beams to RIS
        c = c0;
        c.rng_ind = sim_i;
        c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
        c.fc = 300e9;
        c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
        num_to_ris = 16;
        c.BFmatB(:,1:num_to_ris) = repmat([c.phiBR_loc c.thetaBR_loc]', 1, num_to_ris);
        c.BFmatM(:,1:num_to_ris) = repmat([c.phiMR_loc c.thetaMR_loc]', 1, num_to_ris);
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(2, sim_i) = PEB;
        OEB_cell{iter}(2, sim_i) = OEB;
        
        
        % Realization 3: 8 beams to RIS
%         c = c0;
%         c.rng_ind = sim_i;
%         c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
%         c.fc = 300e9;
%         c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
%         c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
%         c.NsaM_dim = [5 5];    
%         c.DsaB = 5;
%         c.DsaM = 5;
%         c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
%         c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
%         num_to_ris = 8;
%         c.BFmatB(:,1:num_to_ris) = repmat([c.phiBR_loc c.thetaBR_loc]', 1, num_to_ris);
%         c.BFmatM(:,1:num_to_ris) = repmat([c.phiMR_loc c.thetaMR_loc]', 1, num_to_ris);
%         c = updateParameters(c);
%         [FIM] = get_fim_uplink(c);
%         [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
%         PEB_cell{iter}(3, sim_i) = PEB;
%         OEB_cell{iter}(3, sim_i) = OEB;
        
        
        % Realization 4: 0 beams to RIS
        c = c0;
        c.rng_ind = sim_i;
        c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
        c.fc = 300e9;
        c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(4, sim_i) = PEB;
        OEB_cell{iter}(4, sim_i) = OEB;
        
        % Realization 4: quantization 1-bit
        c = c0;
        c.rng_ind = sim_i;
        c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
        c.fc = 300e9;
%         c.RIS_optimizer = "Random";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.RIS_optimizer = "Quantized";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.RIS_quantization = 1;
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
        num_to_ris = 16;
        c.BFmatB(:,1:num_to_ris) = repmat([c.phiBR_loc c.thetaBR_loc]', 1, num_to_ris);
        c.BFmatM(:,1:num_to_ris) = repmat([c.phiMR_loc c.thetaMR_loc]', 1, num_to_ris);
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(5, sim_i) = PEB;
        OEB_cell{iter}(5, sim_i) = OEB;
        
        % Realization 5: quantization 2-bit
        c = c0;
        c.rng_ind = sim_i;
        c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
        c.fc = 300e9;
        c.RIS_optimizer = "Quantized";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.RIS_quantization = 2;
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
        num_to_ris = 16;
        c.BFmatB(:,1:num_to_ris) = repmat([c.phiBR_loc c.thetaBR_loc]', 1, num_to_ris);
        c.BFmatM(:,1:num_to_ris) = repmat([c.phiMR_loc c.thetaMR_loc]', 1, num_to_ris);
        c = updateParameters(c);
        [FIM] = get_fim_uplink(c);
        [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        PEB_cell{iter}(6, sim_i) = PEB;
        OEB_cell{iter}(6, sim_i) = OEB;
        
        % Greedy
        c = c0;
        c.rng_ind = sim_i;
        c.NR_dim = [sqrt(Iter_vec(iter)) sqrt(Iter_vec(iter))]';   
        c.fc = 300e9;
        c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
        c.NsaB_dim = [5 5];     % antenna elements (AEs) per SA at BS
        c.NsaM_dim = [5 5];    
        c.DsaB = 5;
        c.DsaM = 5;
        c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
        c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
        
        PEB_all = zeros(1,16);
        for bf_i = 1:16
            disp(bf_i);
            BS_to_RIS = bf_i;
            MD_to_RIS = bf_i;
            c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
            c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';
            c.BFmatB(:,1:BS_to_RIS) = repmat([c.phiBR_loc c.thetaBR_loc]', 1, BS_to_RIS);
            c.BFmatM(:,1:MD_to_RIS) = repmat([c.phiMR_loc c.thetaMR_loc]', 1, MD_to_RIS);
            c = updateParameters(c);
            [FIM] = get_fim_uplink(c);
            [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
            PEB_all(bf_i) = PEB;
        end
%         figure;plot(PEB_all);
        PEB_cell{iter}(3, sim_i) = min(PEB_all);
        OEB_cell{iter}(3, sim_i) = OEB;

    end
end

% save data_Fig_5.mat
% load data_Fig_5.mat
%% plot
% PEB_data = zeros(6, length(Iter_vec));
% OEB_data = zeros(6, length(Iter_vec));
% 
% for i = 1:length(Iter_vec)
%     PEB_data(:, i) = mean(PEB_cell{i},2);
%     OEB_data(:, i) = mean(OEB_cell{i},2);
% end
Iter_vec = (1:20).^2;
linewidth = 1.5;
markersize = 8;
% PEB
figure;semilogy(Iter_vec, PEB_data(4,:), 'k-', 'Linewidth', linewidth, 'MarkerSize', markersize);
hold on;semilogy(Iter_vec, PEB_data(5,:), 'rs--', 'Linewidth', linewidth, 'MarkerSize', markersize);
hold on;semilogy(Iter_vec, PEB_data(6,:), 'm--', 'Linewidth', linewidth, 'MarkerSize', markersize);
hold on;semilogy(Iter_vec, PEB_data(2,:), 'g+-', 'Linewidth', linewidth, 'MarkerSize', markersize);
% hold on;semilogy(Iter_vec, PEB_data(4,:), 'm^-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, PEB_data(3,:), 'bo-', 'Linewidth', linewidth, 'MarkerSize', markersize);

legend('{THz Benchmark ($b_\mathrm{R}$=$0$)}', ...
       '{THz  ($b_\mathrm{R}$=$16$, QTZ=1 bit)}', ...
       '{THz  ($b_\mathrm{R}$=$16$, QTZ=2 bit)}', ...
       '{THz  ($b_\mathrm{R}$=$16$)}', ...
       '{THz (Adaptive $b_\mathrm{R}$)}',...
       'Location', 'NorthEast', 'Interpreter','Latex');
xlabel('{Number of RIS Elements ($N_\mathrm{R}$)}','Interpreter','Latex')
ylabel('PEB [m]','Interpreter','Latex');

grid on;

set(gca,'fontsize', 18);
set(gcf,'position', [100,100,400*1.3, 400])


% axis([min(Iter_vec) max(Iter_vec) 1e-1 10])
% axis([min(Iter_vec) max(Iter_vec) 0.1 0.6]);

% grid on;
% set(gca,'fontsize', 14);
% set(gcf,'position', [100,100,400*1.3, 400])

% print -dpng -r600 sim_5_RIS_NR_PEB.png

% saveas(gcf,'sim_5_RIS_NR_PEB', 'epsc')

%% OEB [^\circ]
figure;plot(Iter_vec, OEB_data(4,:), 'r--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(5,:), 'g+-', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(3,:), 'ro--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(2,:), 'm+--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(1,:), 'b-', 'Linewidth', 1.5);
hold on;plot(Iter_vec, OEB_data(6,:), 'k-', 'Linewidth', linewidth);


legend('{Random Coefficient}', ...
       '{Fixed Coefficient}', ...
       '{Optimized, 1-bit QNT ($\Delta_{\mathrm{R}}$=$0.5{\lambda})$}', ...
       '{Optimized, 2-bit QNT ($\Delta_{\mathrm{R}}$=$0.5{\lambda})$}', ...
       '{Optimized, w/o QNT ($\Delta_{\mathrm{R}}$=$0.5{\lambda})$}',...
       '{Optimized, w/o QNT ($\Delta_{\mathrm{R}}$=$2.5{\lambda})$}',...
       'Location', 'NorthEast', 'Interpreter','Latex');
   
xlabel('{Number of RIS Elements $N_\mathrm{R}$}','Interpreter','Latex')
ylabel('OEB [deg]','Interpreter','Latex');
% axis([min(Iter_vec) max(Iter_vec) 5 15])

grid on;
set(gca,'fontsize', 14);
set(gcf,'position', [100,100, 400*1.29, 400])

% print -dpng -r600 sim_4_RIS_NR_OEB.png

%%





