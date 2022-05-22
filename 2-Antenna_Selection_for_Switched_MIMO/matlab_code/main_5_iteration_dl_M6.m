close all;
clear all;
clc;
% random doa simulation
%% 
umax = 0.9;
% ref_angle = -thetam:1:thetam;
% ref_doa = sin(ref_angle/180*pi);
ref_doa = -umax:0.01:umax;

thetam = floor(asind(umax));
max_doa = 0.9;
L = 10;
interval = 0.5;
TRA_max = 0.9;
TRA_interval = 0.1;
D = 0:interval:L;
% snrdB = -40
% crlb = calculate_crlb(S1, snr);

sim_times = 5000;
rng(0);

%%
% load sim_M6_0915.mat 
load sim_M6_1124.mat 
load model_continuous_M6_1123.mat;
% load model_discrete_M6_1123.mat;

M = 6;
snrdB = -5:2.5:25; % M6

tic
for i = 1:length(snrdB)
    snr = snrdB(i)

    S7 = alg_nn(abs(sin(target_doa/180*pi)), (snr/30-0.5)*2, variables, D, M);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S7', snr);
    MSE_avg(7,i) = Es(2);
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution

    % NN
    for j = 1:sim_times
        S7 = alg_nn(abs(sin(target_all(j)/180*pi)), (snr/30-0.5)*2, variables, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S7, snr);
    end
%     [results, outputs] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim(7,i) = mean((results).^2);
    
    toc
end
%% Approximation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_avg(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_avg(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_avg(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_avg(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('CRLB-PSL [\beta = 1]', 'CRLB-PSL [\beta = 0.75]', 'CRLB-PSL [\beta = 0.5]', ...
'ULA', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)
%% Simulation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_sim(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_sim(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('PSL-C [\beta = 1]', 'PSL-C [\beta = 0.75]', 'PSL-C [\beta = 0.5]', ...
'ULA', 'TRA-G', 'TRA-DL', 'Best Case CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)

% save sim_data_M6_1124.mat
% load sim_data_M6_1124.mat

% print -dpng -r600 sim-M6.png

