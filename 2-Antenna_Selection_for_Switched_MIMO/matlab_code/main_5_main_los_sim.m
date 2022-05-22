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

sim_times = 50;
rng(0);

%% comparison
M = 8;
snrdB = -5:2.5:25;  % M8

% Benchmark curves: 1-CRLB-PSL-1; 2-0.75; 3-0.5; 4-ULA; 5-TRA-Greedy; 
% 6-TRA-Greedy-p; 7-TRA-NN-p; 8-Optimal CRLB

MSE_avg = zeros(8, length(snrdB));
MSE_sim = zeros(8, length(snrdB));
target_doa = 0;

% benchmark 1 to 4
delta = 1;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S1, ~] = linear_iso_search(D,M,delta);

delta = .75;
% [P_am1, x_am] = linear_amb(D,M,delta);    % not as good as search.
[S2, ~] = linear_iso_search(D,M,delta);

delta = 0.5;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S3, ~] = linear_iso_search(D,M,delta);

S4 = 0:(M-1);



tic
for i = 1:length(snrdB)
    snr = snrdB(i)
    
    % Sol 8: CRLB
    S8 = alg_crlb_greedy(snr, D, M);
    
    % Sol 6: threshold greedy. single point
    S6 = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);
    
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S1', snr);
    MSE_avg(1,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
    MSE_avg(2,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
    MSE_avg(3,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S4', snr);
    MSE_avg(4,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S6', snr);
    MSE_avg(6,i) = Es(2);
    
    
    crlb = calculate_crlb(S8, snr);
    MSE_avg(8,i) = crlb;
    
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim(1,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim(2,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim(3,i) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim(4,i) = mean((results).^2);
    
    
    MSE_sim(8,i) = crlb;


    % TRA-greedy-p
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S6 = alg_threshold_greedy(target_all(j), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S6, snr);
    end
%     [results, outputs] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim(6,i) = mean((results).^2);

    toc
end

save sim_M8.mat 
% load sim_M8.mat 
%% Approximation
% linewidth = 1.2;
% markersize = 8;
% figure;semilogy(snrdB, MSE_avg(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% % hold on;semilogy(snrdB, MSE_avg(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_avg(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
% hold on;semilogy(snrdB, MSE_avg(7,:), '--b>', 'Linewidth', linewidth, 'MarkerSize',markersize)
% 
% % hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
% 
% 
% grid on;
% legend('PSLC [\beta = 1]', 'PSLC [\beta = 0.75]', 'PSLC [\beta = 0.5]', ...
% 'ULA', 'TRA-G', 'TRA-DL', 'Best Case CRLB')
% xlabel('SNR (dB)')
% ylabel('MSE')
% grid on; set(gca,'FontSize',12)
%% Simulation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_sim(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '--oc', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_sim(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('PSLC [\beta = 1]', 'PSLC [\beta = 0.75]', 'PSLC [\beta = 0.5]', ...
'ULA', 'TRA-G', 'TRA-DL', 'Best Case CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)

%% comparison
M = 6;
snrdB = -5:2.5:25; % M6

% Benchmark curves: 1-CRLB-PSL-1; 2-0.75; 3-0.5; 4-ULA; 5-TRA-Greedy; 
% 6-TRA-Greedy-p; 7-TRA-NN-p; 8-Optimal CRLB

MSE_avg = zeros(8, length(snrdB));
MSE_sim = zeros(8, length(snrdB));
target_doa = 0;

% benchmark 1 to 4
delta = 1;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S1, ~] = linear_iso_search(D,M,delta);

delta = .75;
% [P_am1, x_am] = linear_amb(D,M,delta);    % not as good as search.
[S2, ~] = linear_iso_search(D,M,delta);

delta = 0.5;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S3, ~] = linear_iso_search(D,M,delta);

S4 = 0:(M-1);



tic
for i = 1:length(snrdB)
    snr = snrdB(i)
    
     % Sol 8: CRLB
    S8 = alg_crlb_greedy(snr, D, M);
    
    % Sol 6: threshold greedy. single point
    S6 = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);
    
  
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S1', snr);
    MSE_avg(1,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
    MSE_avg(2,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
    MSE_avg(3,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S4', snr);
    MSE_avg(4,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S6', snr);
    MSE_avg(6,i) = Es(2);
    
    
    crlb = calculate_crlb(S8, snr);
    MSE_avg(8,i) = crlb;
    
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim(1,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim(2,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim(3,i) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim(4,i) = mean((results).^2);
    
    MSE_sim(8,i) = crlb;


    % TRA-greedy-p
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S6 = alg_threshold_greedy(target_all(j), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S6, snr);
    end
%     [results, outputs] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim(6,i) = mean((results).^2);

    toc
end

% save simulation_M8.mat 
% load simulation_M8.mat 

% save simulation_M6.mat 
% load simulation_M6.mat 

save sim_M6.mat 
% load sim_M6_0804.mat 
% load sim_M6_0915.mat
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
legend('CRLB-PSL [\beta = 1]', 'CRLB-PSL [\beta = 0.75]', 'CRLB-PSL [\beta = 0.5]', ...
'ULA', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)

%%
M = 4;
snrdB = 0:2.5:30; % M4

% Benchmark curves: 1-CRLB-PSL-1; 2-0.75; 3-0.5; 4-ULA; 5-TRA-Greedy; 
% 6-TRA-Greedy-p; 7-TRA-NN-p; 8-Optimal CRLB

MSE_avg = zeros(8, length(snrdB));
MSE_sim = zeros(8, length(snrdB));
target_doa = 0;

% benchmark 1 to 4
delta = 1;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S1, ~] = linear_iso_search(D,M,delta);

delta = .75;
% [P_am1, x_am] = linear_amb(D,M,delta);    % not as good as search.
[S2, ~] = linear_iso_search(D,M,delta);

delta = 0.5;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S3, ~] = linear_iso_search(D,M,delta);

S4 = 0:(M-1);



tic
for i = 1:length(snrdB)
    snr = snrdB(i)
    
    % Sol 8: CRLB
    S8 = alg_crlb_greedy(snr, D, M);
    
    % Sol 6: threshold greedy. single point
    S6 = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S1', snr);
    MSE_avg(1,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
    MSE_avg(2,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
    MSE_avg(3,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S4', snr);
    MSE_avg(4,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S6', snr);
    MSE_avg(6,i) = Es(2);
    
    crlb = calculate_crlb(S8, snr);
    MSE_avg(8,i) = crlb;
    
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim(1,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim(2,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim(3,i) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim(4,i) = mean((results).^2);
    
    MSE_sim(8,i) = crlb;


    % TRA-greedy-p
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S6 = alg_threshold_greedy(target_all(j), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S6, snr);
    end
%     [results, outputs] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim(6,i) = mean((results).^2);

    toc
end

% save simulation_M8.mat 
% load simulation_M8.mat 

% save simulation_M6.mat 
% load simulation_M6.mat 

save sim_M4.mat 
% load sim_M4.mat 

% load sim_M6_0803.mat 
% load sim_M4_1223.mat
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
'ULA', 'TRA-Greedy', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
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
legend('CRLB-PSL [\beta = 1]', 'CRLB-PSL [\beta = 0.75]', 'CRLB-PSL [\beta = 0.5]', ...
'ULA', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)