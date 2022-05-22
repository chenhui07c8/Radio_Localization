% Code for reproducible results in paper:
% "Channel Model Mismatch Analysis for XL-MIMO Systems from a Localization
% Perspective"
% 
% Version: 08-May
% author: Hui Chen (hui.chen@chalmers.se; hui.chen@kaust.edu.sa)
% 
%%
close all;
clear all;
clc;

%% initialization...
% p: true model (near-field model, perfect knowledge of model)
% f: mismatched model (far-field model, approximated model)

rng(2);
c = default_THz_2D_SWM_parameters();
c.array_structure = "Digital";     % Digital
c.wave_type = "PWM";    % PWM, SWM, SWM_A1
c.K = 10;
c.G = 1;
c.P = 100;

c.D_Rx = (0:63)';  % uniform linear array
c.PU = [2 2]'*1;
c.PU0 = c.PU;
c.AOA0 = atan2d(c.PU0(2), c.PU0(1));
c.d0 = norm(c.PU0);
c = update_parameters(c);

c0 = c;
%% PWM model and SWM model.
% cf is the mismatched model using PWM
cf = c0;
cf.wave_type = "PWM";
% cf.wave_type = "SWM";
cf = update_parameters(cf);
cf = get_Rx_symbol(cf);
cf0 = cf;

% cp is the true model using SWM
cp = c;
cp.wave_type = "SWM";
% spatial non-stationarity, spherical wave model, beam squint effect
cp.NF_SNS = "True";     % True, False
cp.NF_SWM = "True";
cp.NF_BSE = "True";
cp = update_parameters(cp);
cp = get_Rx_symbol(cp);
cp.yp = cp.u;
cp0 = cp;


snr = c.snr;
N_iter = 15;
P_vec = (10.^(linspace(-1, 3, N_iter)));
P_vecdB = pow2db(P_vec);

AEB_F = zeros(1, length(P_vec));
DEB_F = zeros(1, length(P_vec));
PEB_F = zeros(1, length(P_vec));

AEB_P = zeros(1, length(P_vec));
DEB_P = zeros(1, length(P_vec));
PEB_P = zeros(1, length(P_vec));

AEB_LB = zeros(1, length(P_vec));
DEB_LB = zeros(1, length(P_vec));
PEB_LB = zeros(1, length(P_vec));

for i = 1:length(P_vec)
    cf.P = P_vec(i);
    cf = get_Rx_symbol(cf);
    cf = get_CRLB(cf);
    AEB_F(i) = cf.AEB;
    DEB_F(i) = cf.DEB;
    PEB_F(i) = cf.PEB;
    
    cp.P = P_vec(i);
    cp = get_Rx_symbol(cp);
    cp = get_CRLB(cp);
    AEB_P(i) = cp.AEB;
    DEB_P(i) = cp.DEB;
    PEB_P(i) = cp.PEB;
    
    [AEB_LB(i), DEB_LB(i), PEB_LB(i)] = get_CRLB_PWM_LB(cp, cf);

end

figure;semilogy(P_vecdB, AEB_F, 'bo--');
hold on;semilogy(P_vecdB, DEB_F, 'ro--');
hold on;semilogy(P_vecdB, PEB_F, 'go--');
% 
hold on;semilogy(P_vecdB, AEB_P, 'b-');
hold on;semilogy(P_vecdB, DEB_P, 'r-');
hold on;semilogy(P_vecdB, PEB_P, 'g-');

hold on;semilogy(P_vecdB, PEB_LB, 'g-');
% legend('AEB-PWM', 'DEB-PWM', 'PEB-PWM');
legend('AEB-PWM', 'DEB-PWM', 'PEB-PWM','AEB-PWM', 'DEB-PWM', 'PEB-PWM');
% LB
%% grid search and crlb
rng(0);
simtimes = 500;

PU_est_mat = zeros(2, simtimes);
DOA_est_vec = zeros(1, simtimes);
delay_est_vec = zeros(1, simtimes);

results_cell = cell(1, length(P_vec));
tic;
parfor i = 1:length(P_vec)
    % i = 1
    % i = length(P_vec)
    cf = cf0;
    cf.P = P_vec(i);
    cf = get_Rx_symbol(cf);
    y0f = cf.u;
    
    cp = cp0;
    cp.P = P_vec(i);
    cp = get_Rx_symbol(cp);
    y0 = cp.u;

    disp(num2str(i) + " of " + num2str(17));
    for sim_i = 1:simtimes
        
        y = y0 + (randn(size(y0)) + 1j*randn(size(y0)))*sqrt(db2pow(-snr)/2);
        angle_resolution = 0.2;
        delay_resolution = 0.02;
        cf.AOA_candidate = c0.AOA0 + randn(1,1)*0.2 + (-10:angle_resolution:10);
        cf.d_candidate = c0.d0 + randn(1,1)*0.02 + [-1:delay_resolution:1];   % initial position
        initial = get_pos_coarse(y, cf);
        initial = initial(1:2);
        cf.PU = initial;
        % debug with an arbitrary initial position
        % cf.PU = c0.PU0+randn(size(c0.PU0))*0.01; 
        % MMLE
        PU_est = alg_mle_toy_POS(y, cf);
        angle_est = atan2d(PU_est(2), PU_est(1));
        d_est = norm(PU_est);
        results_cell{i}(1, sim_i) = norm(c0.AOA0 - angle_est);
        results_cell{i}(2, sim_i) = norm(c0.d0 - d_est);
        results_cell{i}(3, sim_i) = norm(c0.PU0 - PU_est);

        % MLE
        cp.PU = initial;
        PU_est = alg_mle_toy_POS(y, cp);
        angle_est = atan2d(PU_est(2), PU_est(1));
        d_est = norm(PU_est);
        results_cell{i}(4, sim_i) = norm(c0.AOA0 - angle_est);
        results_cell{i}(5, sim_i) = norm(c0.d0 - d_est);
        results_cell{i}(6, sim_i) = norm(c0.PU0 - PU_est);

    end
end
toc;
%%
% figure;plot(results_cell{i}(1, :))
results_mat = zeros(length(P_vec), 6);
for i = 1:length(P_vec)
    for j = 1:size(results_mat, 2)
        results_mat(i, j) = sqrt(mean(results_cell{i}(j,:).^2));
    end
end

%%
% save Fig-1.mat
% load Fig-1.mat

marker_spacing = 1;
linewidth = 2;
markersize = 10;
figure;semilogy(P_vecdB, results_mat(:,3), 'bo--', 'LineWidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:marker_spacing:length(P_vec));
hold on;semilogy(P_vecdB, results_mat(:,6), 'rs--', 'LineWidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:marker_spacing:length(P_vec));
hold on;semilogy(P_vecdB, PEB_LB, 'k+-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'MarkerIndices', 1:marker_spacing:length(P_vec));
hold on;semilogy(P_vecdB, PEB_F, 'b^-', 'LineWidth', linewidth, 'MarkerSize', 8);
hold on;semilogy(P_vecdB, PEB_P, 'r-', 'LineWidth', linewidth);

axis([-10 30 10^-3.2 1e0]);
xlabel('$P$ [dBm]','interpreter','latex');
ylabel('Position RMSE [m]','interpreter','latex');

grid on;
set(gca,'FontSize', 18);
legend('MMLE', 'MLE', ...
    'LB', 'CRB-FF', 'CRB-NF', 'Location', 'Southwest','interpreter','latex');


% addpath('/Users/huiche/Desktop/Projects/Tikz_toolbox/src');matlab2tikz('Fig-1.tex');
