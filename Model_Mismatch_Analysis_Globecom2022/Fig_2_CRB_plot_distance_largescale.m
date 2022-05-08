% Code for reproducible results in paper:
% "Channel Model Mismatch Analysis for XL-MIMO Systems from a Localization
% Perspective"
% 
% Version: 08-May
% author: Hui Chen (hui.chen@chalmers.se; hui.chen@kaust.edu.sa)
% 
%% close all;
clear all;
clc;

%% initialization...
% p: near-field model (assumed true model)
% f: far-field model (mismatched model)
rng(46);

c = default_THz_2D_SWM_parameters();
c.BW = 100e6;
c.array_structure = "Digital";     % Digital
c.wave_type = "PWM";
c.K = 10;
c.G = 1;
c.P = 100;
c.D_Rx = (0:63)';  % uniform linear array


c = update_parameters(c);
c0 = c;
%% PWM model and SWM model.
% cf is the mismatched model using PWM
cf = c0;
cf.wave_type = "PWM";
cf = update_parameters(cf);
cf = get_Rx_symbol(cf);

% cp is the true model using SWM
cp = c0;
cp.wave_type = "SWM";
cp = update_parameters(cp);
cp = get_Rx_symbol(cp);
cp.yp = cp.u;
% c.u

norm(cf.u(:) - cp.u(:), 'fro')/norm(cp.u(:), 'fro')

%% grid search and crlb
% rng(0);
snr = c.snr;
N_iter = 100;

iter_vec = 10.^linspace(-1, 1, N_iter);

% xgrid = 1:1:50;
ygrid = 1;

CRLB_cell = cell(1, length(iter_vec));

P_vec = (c.P);
eta_mat = zeros(4, N_iter);

for xi = 1:length(iter_vec)
    disp([num2str(xi) '/' num2str(length(iter_vec))])
        pos = iter_vec(xi)*[2 2]'/norm([2 2]);
        % PWM
        cf = c0;
        cf.PU = pos;
        cf.wave_type = "PWM";
        cf = update_parameters(cf);
        cf = get_Rx_symbol(cf);
        cf = get_CRLB(cf);
        CRLB_cell{xi}(1, 1) = cf.PEB;
        CRLB_cell{xi}(1, 2) = cf.PEB;
        CRLB_cell{xi}(1, 3) = cf.PEB;
        CRLB_cell{xi}(1, 4) = cf.PEB;
        
        
        % True model consider: with 3 mismatch
        cp = c0;
        cp.PU = pos;
        cp.wave_type = "SWM";
        cp.NF_SNS = "True";
        cp.NF_SWM = "True";
        cp.NF_BSE = "True";
        cp = update_parameters(cp);
        cp = get_Rx_symbol(cp);
        cp = get_CRLB(cp);
        CRLB_cell{xi}(2, 1) = cp.PEB;
        CRLB_cell{xi}(5, 1) = cp.AEB;
        CRLB_cell{xi}(8, 1) = cp.DEB;
        % LB
        [AEB, DEB, PEB] = get_CRLB_PWM_LB(cp, cf);
        CRLB_cell{xi}(3, 1) = PEB;
        CRLB_cell{xi}(6, 1) = AEB;
        CRLB_cell{xi}(9, 1) = DEB;
        % pow2db(abs(PEB-cp.PEB)/cp.PEB)
        
        % True model consider: ignore 1
        cp = c0;
        cp.PU = pos;
        cp.wave_type = "SWM";
        cp.NF_SNS = "True";
        cp.NF_SWM = "False";
        cp.NF_BSE = "False";
        cp = update_parameters(cp);
        cp = get_Rx_symbol(cp);
        cp = get_CRLB(cp);
        CRLB_cell{xi}(2, 2) = cp.PEB;
        CRLB_cell{xi}(5, 2) = cp.AEB;
        CRLB_cell{xi}(8, 2) = cp.DEB;
        % LB
        [AEB, DEB, PEB, eta_opt] = get_CRLB_PWM_LB(cp, cf);
        eta_mat(:, xi) = eta_opt;
        CRLB_cell{xi}(3, 2) = PEB;
        CRLB_cell{xi}(6, 2) = AEB;
        CRLB_cell{xi}(9, 2) = DEB;
        % figure;plot(eta_mat(2,:))
        
        % True model consider: ignore 2
        cp = c0;
        cp.PU = pos;
        cp.wave_type = "SWM";
        cp.NF_SNS = "False";
        cp.NF_SWM = "True";
        cp.NF_BSE = "False";
        cp = update_parameters(cp);
        cp = get_Rx_symbol(cp);
        cp = get_CRLB(cp);
        CRLB_cell{xi}(2, 3) = cp.PEB;
        CRLB_cell{xi}(5, 3) = cp.AEB;
        CRLB_cell{xi}(8, 3) = cp.DEB;
        % LB
        [AEB, DEB, PEB] = get_CRLB_PWM_LB(cp, cf);
        CRLB_cell{xi}(3, 3) = PEB;
        CRLB_cell{xi}(6, 3) = AEB;
        CRLB_cell{xi}(9, 3) = DEB;
        
        % True model consider: ignore 3
        cp = c0;
        cp.PU = pos;
        cp.wave_type = "SWM";
        cp.NF_SNS = "False";
        cp.NF_SWM = "False";
        cp.NF_BSE = "True";
        cp = update_parameters(cp);
        cp = get_Rx_symbol(cp);
        cp = get_CRLB(cp);
        CRLB_cell{xi}(2, 4) = cp.PEB;
        CRLB_cell{xi}(5, 4) = cp.AEB;
        CRLB_cell{xi}(8, 4) = cp.DEB;
        % LB
        [AEB, DEB, PEB] = get_CRLB_PWM_LB(cp, cf);
        CRLB_cell{xi}(3, 4) = PEB;
        CRLB_cell{xi}(6, 4) = AEB;
        CRLB_cell{xi}(9, 4) = DEB;
end
%%


CRB_mat = zeros(length(iter_vec), 9, 4);

for ti = 1:size(CRB_mat, 2)
    for xi = 1:length(iter_vec)
        CRB_mat(xi, ti, :) = CRLB_cell{xi}(ti, :);
    end
end

size(CRB_mat)

%% PEB
PEB_PWM = reshape(CRB_mat(:, 1, :), [length(iter_vec), 4]);
PEB_SWM = reshape(CRB_mat(:, 2, :), [length(iter_vec), 4]);
PEB_LB = reshape(CRB_mat(:, 3, :), [length(iter_vec), 4]);

% figure;plot(PEB_LB)

E_PEB = pow2db(abs(PEB_LB - PEB_SWM)./PEB_SWM);
% E_PEB = (abs(PEB_LB - PEB_SWM));

figure;semilogx(iter_vec, E_PEB(:,1), 'k-', 'Linewidth', 2);
hold on;plot(iter_vec, E_PEB(:,2), 'bo-', 'MarkerIndices', 1:4:length(iter_vec));
hold on;plot(iter_vec, E_PEB(:,3), 'rs-', 'MarkerIndices', 1:4:length(iter_vec));
hold on;plot(iter_vec, E_PEB(:,4), 'gd-', 'MarkerIndices', 1:4:length(iter_vec));
xlabel('Distance [m]');
ylabel('Relative Mismatch Error [dB]');
legend('LB (TM)', 'LB (only SNS)', 'LB (only SWM)', 'LB (only BSE)', 'Location', 'Northeast');
axis([0.25 10 -60 20]);
grid on;
set(gca,'FontSize', 18);
%% angle
PEB_PWM = reshape(CRB_mat(:, 4, :), [length(iter_vec), 4]);
PEB_SWM = reshape(CRB_mat(:, 5, :), [length(iter_vec), 4]);
PEB_LB = reshape(CRB_mat(:, 6, :), [length(iter_vec), 4]);

E_PEB = pow2db(abs(PEB_LB - PEB_SWM)./PEB_SWM);

figure;semilogx(iter_vec, E_PEB(:,1), 'k-', 'Linewidth', 2);
hold on;plot(iter_vec, E_PEB(:,2), 'bo--', 'MarkerIndices', 1:4:length(iter_vec));
hold on;plot(iter_vec, E_PEB(:,3), 'rs--', 'MarkerIndices', 1:4:length(iter_vec));
hold on;plot(iter_vec, E_PEB(:,4), 'gd--', 'MarkerIndices', 1:4:length(iter_vec));
xlabel('Distance [m]');
ylabel('Relative Mismatch Error [dB]');
% legend('LB (TM)', 'LB (only SNS)', 'LB (only SWM)', 'LB (only BSE)', 'Location', 'Northeast');
% axis([0 150 -80 0]);
axis([0.25 10 -60 30]);
grid on;
set(gca,'FontSize', 18);



%% delay
PEB_PWM = reshape(CRB_mat(:, 7, :), [length(iter_vec), 4]);
PEB_SWM = reshape(CRB_mat(:, 8, :), [length(iter_vec), 4]);
PEB_LB = reshape(CRB_mat(:, 9, :), [length(iter_vec), 4]);

E_PEB = pow2db(abs(PEB_LB - PEB_SWM)./PEB_SWM);

figure;semilogx(iter_vec, E_PEB(:,1), 'k-', 'Linewidth', 2);
hold on;plot(iter_vec, E_PEB(:,2), 'bo--', 'MarkerIndices', 1:4:length(iter_vec));
hold on;plot(iter_vec, E_PEB(:,3), 'rs--', 'MarkerIndices', 1:4:length(iter_vec));
hold on;plot(iter_vec, E_PEB(:,4), 'gd--', 'MarkerIndices', 1:4:length(iter_vec));
xlabel('Distance [m]');
ylabel('Relative Mismatch Error [dB]');
% legend('LB (TM)', 'LB (only SNS)', 'LB (only SWM)', 'LB (only BSE)', 'Location', 'Northeast');
% axis([0 150 -80 0]);
axis([0.25 10 -60 20]);

grid on;
set(gca,'FontSize', 18);

