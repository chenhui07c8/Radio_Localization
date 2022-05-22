close all;
clear all;
clc;

%% Initialization
c = defaultParameters();
c.rng_ind = 0;
% BM & BRM channel
c.NB_dim = [1 20]';   % Dim of BS Elements
c.NR_dim = [4 4]';   
c.NM_dim = [4 4]';
c.PB = [0 0 0]';
c.PM = [10, 0, 0]';
c.OB = [0, 0, 0]';
c.OR = [-90, 0, 0]';
c.OM = [180, 0, 0]';

c.NM_dim = [1 1]';
c.PC = [];  % assign an empty matrix if not available
c.PR = [];  % assign an empty matrix if not available

% signal parameters
% mmWave
% c.fc = 60e9;
% c.BW = 100e6;    % 100M bandwidth
% c.K = 10;

% THz
c.fc = 300e9;
c.BW = 100e6;    % 100M bandwidth
c.K = 10;


% Hardware Parameters
% c.P = 10;   % Transmission power: 1--> 1 dBm
% c.GB = [1 180 180];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;
% c.GM = [1 180 180];  % [gain, max_phi, max_theta], HPBW = 2*max_phi;

% SA: analog beamforming
c.deltaDaeCoe = 1;    % in [lambda/2]. coe=1 --> 0.5*lambdac.

sa_dim_B = 5;     % antenna elements (AEs) per SA at BS
sa_dim_M = 1;
c.NsaB_dim = [sa_dim_B sa_dim_B];     % antenna elements (AEs) per SA at BS
c.NsaM_dim = [sa_dim_M sa_dim_M];
c.BFmatB = zeros(2, prod(c.NB_dim))+[0 0]';  % row1: phi, row2: theta
c.BFmatM = zeros(2, prod(c.NM_dim))+[0 0]';
c.DsaB = sa_dim_B;
c.DsaM = sa_dim_M;



% Model parameters
c.symbol_type = "Random"; % Ones, Random, RandomSymbol
c.pos_type = "2D";         % 2D, 3D positions
c.ori_type = "0D";
c.syn_type = "Asyn"; % Syn, Asyn

% c.AoSA_type = "SSP"; % SWM/PWM model for: SA phase/SA amplitude/AE
c.AoSA_type = "SPP"; % SWM/PWM model for: SA phase/SA amplitude/AE

c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
% c.RIS_optimizer = "Random";   % Elzanaty, He, Random, Quantized, Fixed, Customized
% c.fixed_noise = "true";       % false, true
% c.operationBW = 10e6;

% different realizations
c = updateParameters(c);


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
c = c0;
c.fc = 60e9;
c.BW = 100e6;    % 100M bandwidth
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));

c = c0;
c.fc = 300e9;
c.BW = 1e8;    % 100M bandwidth
c = updateParameters(c);
[FIM] = get_fim_uplink(c);
[CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
disp("PEB = " + num2str(PEB) + " | " + "OEB = " + num2str(OEB));

%% SIM-1: SA size (lambda/2, lambda/8)
% close all;
step = 20;
Iter_vec = 10.^linspace(-1,1,step);
PEB_cell = cell(1, length(Iter_vec));
PEB_PWM_cell = cell(1, length(Iter_vec));
DEB_cell = cell(1, length(Iter_vec));
AEB_cell = cell(1, length(Iter_vec));


parfor iter = 1:length(Iter_vec)
    disp("Iteration: " + num2str(iter));

    % Realization 1: mmWave
    c = c0;
%     c.BW = Iter_vec(iter)*1e9;
    c.PM = [Iter_vec(iter) 0 0]';
    c.syn_type = "Asyn"; % Syn, Asyn
    c.channel_knowledge = "False";  % True , False, Partial
    c.K = 10;
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(1) = PEB;
    
    [FIM] = get_FIM_PWM(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_PWM_cell{iter}(1) = PEB;
    
    
    [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c);
    DEB_cell{iter}(1) = DEB;
    AEB_cell{iter}(1) = AEB;
    
    % Realization 2: THz, Syn
    c = c0;
%     c.BW = Iter_vec(iter)*1e9;
    c.PM = [Iter_vec(iter) 0 0]';
    c.syn_type = "Asyn"; % Syn, Asyn
    c.channel_knowledge = "Partial";  % True , False, Partial
    c.K = 10;
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(2) = PEB;
    
    [FIM] = get_FIM_PWM(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_PWM_cell{iter}(2) = PEB;
    
    [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c);
    DEB_cell{iter}(2) = DEB;
    AEB_cell{iter}(2) = AEB;

    
    % Realization 3: THz, Syn, AOSA
    c = c0;
    c.PM = [Iter_vec(iter) 0 0]';
    c.syn_type = "Asyn"; % Syn, Asyn
    c.channel_knowledge = "True";  % True , False, Partial
    c.K = 10;
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(3) = PEB;

    [FIM] = get_FIM_PWM(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_PWM_cell{iter}(3) = PEB;
    
    [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c);
    DEB_cell{iter}(3) = DEB;
    AEB_cell{iter}(3) = AEB;



   % Realization 1: mmWave
    c = c0;
    c.PM = [Iter_vec(iter) 0 0]';
    c.syn_type = "Syn"; % Syn, Asyn
    c.channel_knowledge = "False";  % True , False, Partial
    c.K = 10;
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(4) = PEB;
    
    [FIM] = get_FIM_PWM(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_PWM_cell{iter}(4) = PEB;
    
    [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c);
    DEB_cell{iter}(4) = DEB;
    AEB_cell{iter}(4) = AEB;
    
    
    % Realization 2: THz, Syn
    c = c0;
    c.PM = [Iter_vec(iter) 0 0]';
    c.syn_type = "Syn"; % Syn, Asyn
    c.channel_knowledge = "Partial";  % True , False, Partial
    c.K = 10;
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(5) = PEB;
    
    [FIM] = get_FIM_PWM(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_PWM_cell{iter}(5) = PEB;  
    
    [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c);
    DEB_cell{iter}(5) = DEB;
    AEB_cell{iter}(5) = AEB; 
    
    % Realization 3: THz, Syn, AOSA
    c = c0;
    c.PM = [Iter_vec(iter) 0 0]';
    c.syn_type = "Syn"; % Syn, Asyn
    c.channel_knowledge = "True";  % True , False, Partial
    c.K = 10;
    c = updateParameters(c);
    [FIM] = get_fim_uplink(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_cell{iter}(6) = PEB;
    
    [FIM] = get_FIM_PWM(c);
    [CRLB, PEB, OEB] = get_crlb_from_fim(FIM, c);
    PEB_PWM_cell{iter}(6) = PEB;
    
    [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c);
    DEB_cell{iter}(6) = DEB;
    AEB_cell{iter}(6) = AEB;
    
end

%% Plot
% save data_Fig_3.mat

% load data_Fig_3.mat

PEB_data = zeros(6, length(Iter_vec));
PEB_PWM_data = zeros(6, length(Iter_vec));
DEB_data = zeros(6, length(Iter_vec));
AEB_data = zeros(6, length(Iter_vec));

for i = 1:length(Iter_vec)
    PEB_data(:, i) = PEB_cell{i}';
    PEB_PWM_data(:, i) = PEB_PWM_cell{i}';
    DEB_data(:, i) = DEB_cell{i}';
    AEB_data(:, i) = AEB_cell{i}';
end


linewidth = 1.7;
markersize = 12;

markerindice = 1;

% PEB_data = abs(PEB_data);
% PEB

figure;loglog(Iter_vec, PEB_data(1,:), 'k-', 'Linewidth', 1.5, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_data(2,:), 'g+-', 'Linewidth', 2, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_data(3,:), 'rs-', 'Linewidth', linewidth, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_data(4,:), 'r>-', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_data(5,:), 'bo-', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_data(6,:), 'm-', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));

% hold on;loglog(Iter_vec, PEB_PWM_data(1,:), 'k--', 'Linewidth', 1.5);
hold on;loglog(Iter_vec, PEB_PWM_data(2,:), 'g+-.', 'Linewidth', 2, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_PWM_data(3,:), 'rs--', 'Linewidth', linewidth, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_PWM_data(4,:), 'r>--', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_PWM_data(5,:), 'bo--', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_PWM_data(6,:), 'go--', 'Linewidth', linewidth);

grid on;
axis([1e-1 10 5e-6 5]);
xlabel('{Distance [m]}','Interpreter','Latex')
ylabel('PEB [m]','Interpreter','Latex');

set(gca,'fontsize', 18);

% print -dpng -r600 sim_3_SWM_PWM_PEB.png

% Zoomed plot ****************************************
bounds = [0.8 1 5*10^-3.5 5*10^-2.8]; % area to be zoomed
pos = [0.23 0.65 0.26 0.2];      % projected area [x, y, x_width, y_width]
vertex = [2 4];                 % vertex to be linked
p = gca;

% Calculate x,y points of zoomPlot
x1 = 10^((pos(1)-p.Position(1))/p.Position(3)*diff(log10(p.XLim))+log10(p.XLim(1)));
x2 = 10^((pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(log10(p.XLim))+log10(p.XLim(1)));
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));

rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
offset = [0 0 0 0];
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2 + offset(1)], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2 + offset(2)], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1 + offset(3)], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1 + offset(4)], '-.k'); % Line to vertex 3
end
% legend hide
z = axes('position',pos);

box on % put box around new pair of axes
loglog(Iter_vec, PEB_data(1,:), 'k-', 'Linewidth', 1.5, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_data(2,:), 'g+-', 'Linewidth', 2, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_data(3,:), 'rs-', 'Linewidth', linewidth, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_data(4,:), 'r>-', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_data(5,:), 'bo-', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_data(6,:), 'm-', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));

% hold on;loglog(Iter_vec, PEB_PWM_data(1,:), 'k--', 'Linewidth', 1.5);
hold on;loglog(Iter_vec, PEB_PWM_data(2,:), 'g+-.', 'Linewidth', 2, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_PWM_data(3,:), 'rs--', 'Linewidth', linewidth, 'MarkerSize', 10, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_PWM_data(4,:), 'r>--', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
hold on;loglog(Iter_vec, PEB_PWM_data(5,:), 'bo--', 'Linewidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:markerindice:length(Iter_vec));
% hold on;loglog(Iter_vec, PEB_PWM_data(6,:), 'go--', 'Linewidth', linewidth);

axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})


legend('{Asyn/Unknown/SWM}', '{Asyn/Partial/SWM}', ...
       '{Syn/Unknown/SWM}', '{Syn/Partial/SWM}', ...
       '{Asyn/Partial/PWM}', ...
       '{Syn/Unknown/PWM}', '{Syn/Partial/PWM}', ...
       'Location', 'SouthEast', 'Interpreter','Latex', 'NumColumns', 1, 'Location', 'SouthEast');
set(gca,'fontsize', 16);

set(gcf,'position', [100,100,480*1.4, 480])

dim = [.145 .13 .05 .09];
annotation('ellipse',dim);
x = [0.23 0.19];
y = [0.18 0.18];
annotation('textarrow',x,y,'String','SWM', 'FontSize', 15, 'Interpreter','Latex');

dim = [.145 .24 .05 .23];
annotation('ellipse', dim);
x = [0.22 0.18];
y = [0.54 0.45];
annotation('textarrow',x,y,'String','PWM', 'FontSize', 15, 'Interpreter','Latex')

dim = [.84 .61 .05 .09];    % x, y, circ width, circ height
annotation('ellipse',dim);
x = [0.8 0.85];
y = [0.55 0.62];
annotation('textarrow',x,y,'String','Syn', 'FontSize', 15, 'Interpreter','Latex')

dim = [.84 .72 .05 .20];
annotation('ellipse', dim);
x = [0.8 0.85];
y = [0.85 0.85];
annotation('textarrow',x,y,'String','Asyn', 'FontSize', 15, 'Interpreter','Latex')

% print -dpng -r600 sim_3_SWM_PWM_PEB.png

% saveas(gcf,'sim_3_SWM_PWM_PEB', 'epsc')



