% Matlab Code for Paper Published in IEEE Communication Surveys & Tutorials
% "A Tutorial on Terahertz-Band Localization for 6G Communication Systems"
% Author: Hui Chen
% Email: hui.chen@chalmers.se; hui.chen@kaust.edu.sa

close all;
clear all;
clc;
folder = 'THz Loc CRLB_v1.0';
addpath(genpath(folder));

% rng(0);
% sequence naming priority: k, b, m, r
%% Initialization
c = defaultParameters();

c.rng_ind = 3;
c.NB_dim = [4 4]';   % Dim of BS Elements
c.NR_dim = [4 4]';   
c.NM_dim = [2 2]';

c.PB = [0 0 0]';
c.PR = [2.5, 2.5, 0]';
c.PM = [5, 0, 0]';
c.OB = [0, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [180, 0, 0]';
c.K = 10;
c.G = 1;


% NLOS, reflectors
c.PC = [2.5 -2.5 0]';
c.coeClu = [0.9];



% mmwave
% c.fc = 60e9;    % carrier frequency
% c.BW = 100e6;    % 100M bandwidth
% sa_dim_R = 5;
% c.NsaR_dim = [sa_dim_R sa_dim_R];     % antenna elements (AEs) per SA at BS
% c.DsaR = sa_dim_R;

% scene 1       save data_Fig_7_01.mat
% c.OM = [180, 0, 0]';
% beamsToRIS = 0;
% prior_info = false;

      
% THz
c.fc = 300e9;    % carrier frequency
c.BW = 0.1e9;    % 100M bandwidth
c.K = 10;
sa_dim_B = 5;
sa_dim_M = 5;
c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
c.NsaM_dim = [sa_dim_M sa_dim_M];
c.DsaB = sa_dim_B;
c.DsaM = sa_dim_M;
sa_dim_R = 25;
c.NsaR_dim = [sa_dim_R sa_dim_R];     % antenna elements (AEs) per SA at BS
c.DsaR = sa_dim_R;

% scene 2      save data_Fig_7_02.mat
c.OM = [180, 0, 0]';
beamsToRIS = 0;
prior_info = false;

% scene 3      save data_Fig_7_03.mat
% c.OM = [180, 0, 0]';
% beamsToRIS = 1;
% prior_info = true;

% scene 4      save data_Fig_7_04.mat
% c.OM = [150, 0, 0]';
% beamsToRIS = 1;
% prior_info = true;



% Hardware Parameters
% Antenna gain
c.symbol_type = "Random"; % Ones, Random, RandomSymbol, Customized
c.pos_type = "2D";         % 2D, 3D positions
c.ori_type = "1D";
c.syn_type = "Asyn"; % Syn, Asyn
c.AoSA_type = "SPP"; % SWM/PWM model for: SA phase/SA amplitude/AE
c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized


c = updateParameters(c, "All");


c = updateParameters(c, "All");
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
x_vec = -1:0.1:6;
y_vec = -3.5:0.1:3.5;

% x_vec = -1:0.5:6;
% y_vec = -3.5:0.5:3.5;

sim_times = 1;
PEB_cell = cell(1, length(x_vec));
OEB_cell = cell(1, length(x_vec));

c0.multipleTransmission = "Random";   % Random, Fixed (RIS coefficients)

tic;
parfor x_i = 1:length(x_vec)
    clc;
    disp("Iteration: " + num2str(x_i));
    PEB_cell{x_i}(1:6, 1:sim_times) = 0;
    OEB_cell{x_i}(1:6, 1:sim_times) = 0;
    
    for y_i = 1:length(y_vec)
        disp("Iteration: " + num2str(x_i) + " | Iter y: " + num2str(y_i));
        % Realization 1: benchmark without NLOS
        c = c0;
%         c.rng_ind = 1;
        c.PM = [x_vec(x_i), y_vec(y_i), 0]';
        c.beam_split = "True";          % True, False
        c = updateParameters(c);
        
        % with prior
        if prior_info
            c.BFmatB = zeros(2, prod(c.NB_dim))+[c.phiBM_loc c.thetaBM_loc]';  % row1: phi, row2: theta
            c.BFmatM = zeros(2, prod(c.NM_dim))+[c.phiMB_loc c.thetaMB_loc]';  % row1: phi, row2: theta
            if(c.LR > 0 && beamsToRIS > 0)
                c.BFmatB(:, 1:beamsToRIS*4) = repmat([c.phiBR_loc c.thetaBR_loc]', [1 beamsToRIS*4]);
                c.BFmatM(:, 1:beamsToRIS) = repmat([c.phiMR_loc c.thetaMR_loc]', [1 beamsToRIS]);
                c.BFmatRB = zeros(2,prod(c.NR_dim)) + [c.phiRB_loc c.thetaRB_loc]';
                c.BFmatRM = zeros(2,prod(c.NR_dim)) + [c.phiRM_loc c.thetaRM_loc]';
            end
            c.BFmatB = c.BFmatB + randn(size(c.BFmatB))*2;
            c.BFmatM = c.BFmatM + randn(size(c.BFmatM))*2;
            c.BFmatRB = c.BFmatRB + randn(size(c.BFmatRB))*2;
            c.BFmatRM = c.BFmatRM + randn(size(c.BFmatRM))*2;

            c = updateParameters(c);
            [FIM] = get_fim_uplink(c);
            [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        % random beams
        elseif prior_info == false
            c.multipleTransmission = 'RRR';
            cM = get_multiple_transmission_parameters(c);
            [FIM, ~] = get_FIM_multi_transmission(c, cM);
            [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
        end
        
        PEB_cell{x_i}(1,y_i) = PEB;
        OEB_cell{x_i}(1,y_i) = OEB;
    end
end
toc;


%%
% load data_Fig_7_01.mat

% load data_Fig_7_02.mat

% load data_Fig_7_03.mat

% load data_Fig_7_04.mat

% close all
PEB_data = zeros(length(x_vec), length(y_vec));
OEB_data = zeros(length(x_vec), length(y_vec));

for i = 1:length(x_vec)
    PEB_data(i, :) = PEB_cell{i}(1,:);
    OEB_data(i, :) = OEB_cell{i}(1,:);
end


figure;imagesc(x_vec, y_vec, abs(PEB_data'), 'Interpolation', 'bilinear');
hold on;contour(x_vec, y_vec, abs(PEB_data'), [1e-3 0.01 0.1 1], 'ShowText','on', 'LineColor', 'w');

hold on; plot(c.PB(1), c.PB(2), 'rx', 'Linewidth', 2, 'MarkerSize', 14);
if(c.LR > 0)
    hold on; plot(c.PR(1), c.PR(2), 'rs', 'Linewidth', 2, 'MarkerSize', 8);
end
if(c.LC > 0)
    hold on; plot(c.PC(1,:), c.PC(2,:), 'ro', 'Linewidth', 2, 'MarkerSize', 8);
end
colorbar;
set(gca,'YDir','normal')
xlabel('{x-axis [m]}','Interpreter','Latex')
ylabel('{y-axis [m]}','Interpreter','Latex');
hcb=colorbar;
hcb.Title.String = "[m]";

set(gca,'colorscale','log')
caxis([0.001 1]);
% caxis([0.01 1]);

set(gca,'fontsize', 15);
set(gcf,'position', [100,100,400*1.2, 400]);


%%
% print -dpng -r600 sim_7_1.png

% print -dpng -r600 sim_7_2.png

% print -dpng -r600 sim_7_3.png

% print -dpng -r600 sim_7_4.png

% saveas(gcf,'sim_7_1', 'epsc')

% saveas(gcf,'sim_7_2', 'epsc')

% saveas(gcf,'sim_7_3', 'epsc')

% saveas(gcf,'sim_7_4', 'epsc')

