close all;
clear all;
clc;
% rng(0);
% sequence naming priority: k, b, m, r
%% Initialization
c = defaultParameters();

c.rng_ind = 0;

c.NB_dim = [4 4]';   % Dim of BS Elements
c.NR_dim = [20 20]';  
c.NM_dim = [4 4]';
c.PB = [0 0 0]';
c.PR = [5, 5, 1]';
c.PM = [10, 0, -1]';
% c.PM = [5, 0, -1]';
c.OB = [0, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [150, 0, 0]';
% c.P = 10;
c.K = 10;


% cluster
c.PC = 2*[4 -4 0; 2 7 0; 4.1 -4.1 0]';
c.coeClu = [0.2 0.2 0.5];
% c.PC = [5 -3.000 0]';
% c.coeClu = 2*[0.2];
% c.PC = [];
c.PR = [];

% mmwave
c.fc = 60e9;    % carrier frequency
c.BW = 100e6;    % 100M bandwidth

% THz
% c.fc = 300e9;    % carrier frequency
% c.BW = 100e6;    % 100M bandwidth

sa_dim_B = 1;
sa_dim_M = 1;
% sa_dim_B = 5;
% sa_dim_M = 5;
c.NaeB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
c.NaeM_dim = [sa_dim_M sa_dim_M];
c.IntB = sa_dim_B;
c.IntM = sa_dim_M;
c.deltaDaeCoe = 1;    % in [lambda/2]. coe=1 --> 0.5*lambdac.

% Hardware Parameters
% Antenna gain
% c.GB = [1 pi/2 pi/2];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;
% % c.GR = [1 pi/3 pi/3];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;
% c.GM = [1 pi/2 pi/2];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;

c.deltaDsa = c.lambdac/2;   % assume to be equal for all the SAs (for now)
c.BFmatB = zeros(2, prod(c.NB_dim))+[0 0]';  % row1: phi, row2: theta
c.BFmatM = zeros(2, prod(c.NM_dim))+[0 0]';
% reshape([1 2 3; 4 5 6; 7 8 9]', [1, 9])


% Model Parameters
c.symbol_type = "Random"; % Ones, Random, RandomSymbol
c.pos_type = "2D";         % 2D, 3D positions
c.ori_type = "1D";
c.syn_type = "Asyn"; % Syn, Asyn
c.random_seed = "Fixed";


c.AoSA_type = "SSP"; % SWM/PWM model for: SA phase/SA amplitude/AE
% c.AoSA_type = "SPP"; % SWM/PWM model for: SA phase/SA amplitude/AE

c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
% c.RIS_optimizer = "Random";   % Elzanaty, He, Random, Quantized, Fixed, Customized

    
c = updateParameters(c);

size(c.Xmk)

c0 = c;

%% Simulation
% close all;
Iter_vec = [0.01 0.1:0.1:1];

PEB_cell = cell(6, length(Iter_vec));
OEB_cell = cell(6, length(Iter_vec));
CPEB_cell = cell(6, length(Iter_vec));


l = [5 -5 0; 1 4 0; 9 -4 0; 5.1 -5 0];

parfor iter = 1:length(Iter_vec)
    disp("Iteration: " + num2str(iter));
    PEB_cell{iter}(1:5) = 0;
    OEB_cell{iter}(1:5) = 0;
    CPEB_cell{iter}(1:5) = 0;
    
    
    % Realization 1: benchmark without NLOS
    c = c0;
    c.PC = [];
    c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB, CPEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(1) = PEB;
    OEB_cell{iter}(1) = OEB;
    CPEB_cell{iter}(1) = 1/0;

    % Realization 2: THz, Syn
    c = c0;
%     c.rng_ind = sim_i;
    c.PC = l(1,:)';
    c.coeClu = Iter_vec(iter)*[1];
    c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB, CPEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(2) = PEB;
    OEB_cell{iter}(2) = OEB;
    CPEB_cell{iter}(2) = CPEB(1);
    
    % Realization 2: THz, Syn
    c = c0;
%     c.rng_ind = sim_i;
    c.PC = l([1 2],:)';
    c.coeClu = Iter_vec(iter)*[1 1];
    c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB, CPEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(3) = PEB;
    OEB_cell{iter}(3) = OEB;
    CPEB_cell{iter}(3) = CPEB(1);
   
    % Realization 3: THz, Syn, AOSA
    c = c0;
%     c.rng_ind = sim_i;
    c.PC = l([1:3], :)';
    c.coeClu = Iter_vec(iter)*[1 1 1];
    c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB, CPEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(4) = PEB;
    OEB_cell{iter}(4) = OEB;
    CPEB_cell{iter}(4) = CPEB(1);
    CPEB_cell{iter}(5) = CPEB(2);

    % Realization 4: THz, Syn, AOSA, BF
    c = c0;
%     c.rng_ind = sim_i;
    c.PC = l([1 2 4], :)';
    c.coeClu = Iter_vec(iter)*[1 1 1];
    c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB, CPEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(6) = PEB;
    OEB_cell{iter}(6) = OEB;
    CPEB_cell{iter}(6) = CPEB(3);
   
end


% save data_Fig_6.mat;

% load data_Fig_6.mat;


%% plot
PEB_data = zeros(6, length(Iter_vec));
OEB_data = zeros(6, length(Iter_vec));
CPEB_data = zeros(6, length(Iter_vec));

for i = 1:length(Iter_vec)
    PEB_data(:, i) = PEB_cell{i}';
    OEB_data(:, i) = OEB_cell{i}';
    CPEB_data(:, i) = CPEB_cell{i}';
end


linewidth = 1.2;
% PEB
figure;plot(Iter_vec, PEB_data(1,:), 'k-', 'Linewidth', 1.5);
hold on;plot(Iter_vec, PEB_data(2,:), 'bo--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, PEB_data(3,:), 'rs--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, PEB_data(4,:), 'c+-', 'Linewidth', linewidth);
hold on;plot(Iter_vec, PEB_data(6,:), 'm^-', 'Linewidth', linewidth);

legend('{Only LOS}', '{$\mathcal{L} = \{l_1\}$}', ...
       '{$\mathcal{L} = \{l_1, l_2\}$}',...
       '{$\mathcal{L} = \{l_1, l_2, l_3\}$}',...
       '{$\mathcal{L} = \{l_1, l_2, l_4\}$}',...
       'Location', 'SouthWest', 'Interpreter','Latex');

xlabel('{NLOS Coefficients}','Interpreter','Latex')
ylabel('PEB [m]','Interpreter','Latex');

% axis([0 1 0.03 0.038]);
grid on;
set(gca,'fontsize', 18);
set(gcf,'position', [100,100,350*1.33, 350])


% OEB
figure;plot(Iter_vec, OEB_data(1,:), 'k-', 'Linewidth', 1.5);
hold on;plot(Iter_vec, OEB_data(2,:), 'bo--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(3,:), 'rs--', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(4,:), 'c+-', 'Linewidth', linewidth);
hold on;plot(Iter_vec, OEB_data(6,:), 'm^-', 'Linewidth', linewidth);

legend('{Only LOS}', '{$\mathcal{L} = \{l_1\}$}', ...
       '{$\mathcal{L} = \{l_1, l_2\}$}',...
       '{$\mathcal{L} = \{l_1, l_2, l_3\}$}',...
       '{$\mathcal{L} = \{l_1, l_2, l_4\}$}',...
       'Location', 'SouthWest', 'Interpreter','Latex');

xlabel('{NLOS Coefficients}','Interpreter','Latex')
ylabel('OEB [deg]','Interpreter','Latex');
% axis([0.1 1 0.45 0.65])
grid on;
set(gca,'fontsize', 18);
set(gcf,'position', [100,100,350*1.3, 350])



% RPEB
figure;semilogy(Iter_vec, CPEB_data(6,:), 'm^-', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, CPEB_data(4,:), 'rs--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, CPEB_data(3,:), 'bo--', 'Linewidth', linewidth);
hold on;semilogy(Iter_vec, CPEB_data(2,:), 'k-', 'Linewidth', 1.5);

legend('{$l_1$ $\mathcal{L} = \{l_1, l_2, l_4\}$}', ...
        '{$l_1$ ($\mathcal{L} = \{l_1, l_2, l_3\}$)}',...
        '{$l_1$ ($\mathcal{L} = \{l_1, l_2\}$)}', ...
        '{$l_1$ ($\mathcal{L} = \{l_1\}$)}', ...
       'Location', 'NorthEast', 'Interpreter','Latex');

xlabel('{NLOS Coefficients}','Interpreter','Latex')
ylabel('RPEB [m]','Interpreter','Latex');

axis([0 1 0.05 50])
grid on;
set(gca,'fontsize', 18);
set(gcf,'position', [100,100,350*1.3, 350])



% print -dpng -r600 sim_6_NLOS_PEB.png

% print -dpng -r600 sim_6_NLOS_OEB.png

% print -dpng -r600 sim_6_NLOS_CPEB.png

% print -dpng -r600 sim_6_NLOS_Layout.png


% saveas(gcf,'sim_6_NLOS_PEB', 'epsc')
% saveas(gcf,'sim_6_NLOS_OEB', 'epsc')
% saveas(gcf,'sim_6_NLOS_RPEB', 'epsc')
% saveas(gcf,'sim_6_NLOS_layout', 'epsc')


%%
% close all;
linewidth = 2;
markersize = 8;
figure; plot(c.PB(1), c.PB(2), 'rx',  'Linewidth', linewidth, 'MarkerSize', 12);
hold on; plot(c.PM(1), c.PM(2), 'r^',   'Linewidth', linewidth, 'MarkerSize', markersize);
hold on; plot(l(:,1), l(:,2),'ro',   'Linewidth', linewidth, 'MarkerSize', markersize);
axis([-1 11 -6 6])
hold on; rectangle('Position', [0, -5, 10, 10], 'EdgeColor', 'k', 'LineWidth', 1);

fontsize = 22;
text(c.PB(1) + 0.4, c.PB(2),'BS', 'FontSize', fontsize)
text(c.PM(1) - 1.4, c.PM(2),'UE', 'FontSize', fontsize)
text(l(1,1) - 1, l(1,2)+0.4,'{$l_1$}', 'FontSize', fontsize, 'Interpreter','Latex')
text(l(2,1) +0.4, l(2,2),'{$l_2$}', 'FontSize', fontsize, 'Interpreter','Latex')
text(l(3,1) +0.4, l(3,2) + 0.4,'{$l_3$}', 'FontSize', fontsize, 'Interpreter','Latex')
text(l(4,1) +0.4, l(4,2) - 0.4,'{$l_4$}', 'FontSize', fontsize, 'Interpreter','Latex')

xlabel('{x-axis [m]}','Interpreter','Latex')
ylabel('y-axis [m]','Interpreter','Latex');
legend('BS', 'UE', 'Reflectors','Interpreter','Latex');
set(gca,'fontsize', fontsize);
set(gcf,'position', [100,100,350*1.3, 350])



