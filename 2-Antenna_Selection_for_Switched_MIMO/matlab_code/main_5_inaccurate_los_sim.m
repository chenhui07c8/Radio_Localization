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
L = 10;
interval = 0.5;
TRA_max = 0.9;
D = 0:interval:L;
% snrdB = -40
% crlb = calculate_crlb(S1, snr);
target_doa = 0; % simulation DOA center
max_doa = 0.9;  % simulation DOA max range, sind(65) = 0.9063
tilde_u = 0.05;   % DOA estimation error, 50~75 at 60 deg, -6~6 at 0 deg
% anchors = [-(tilde_u) 0 (tilde_u)]; % 3 anchors to add hat_u
anchors = [-(tilde_u) -(tilde_u)/2 0 (tilde_u)/2 (tilde_u)]; % 3 anchors to add hat_u

sim_times = 5000;
rng(1);

%% comparison
M = 6;
snrdB = -5:2.5:25;  % M8

% Benchmark curves: 1-CRLB-PSL-1; 2-0.85; 3-0.5; 4-ULA; 
% 5-TRA-Greedy []; 6-TRA-Greedy-2; 7-TRA-Greedy-3; 8-Optimal CRLB

MSE_avg_cell = cell(1,length(snrdB));
MSE_sim_cell = cell(1,length(snrdB));
S_all = cell(1,length(snrdB));
for i = 1:length(MSE_avg_cell)
    MSE_avg_cell{i} = zeros(1, 9);
    MSE_sim_cell{i} = zeros(1, 9);
    S_all{i} = zeros(9, M);
end


% benchmark 1 to 4
delta = 1;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S1, ~] = linear_iso_search(D,M,delta);

delta = .85;
% [P_am1, x_am] = linear_amb(D,M,delta);    % not as good as search.
[S2, ~] = linear_iso_search(D,M,delta);

delta = 0.5;
% [P_am2, x_am] = linear_amb(D,M,delta);
[S3, ~] = linear_iso_search(D,M,delta);

S4 = 0:(M-1);



tic
parfor i = 1:length(snrdB)
    snr = snrdB(i);
    disp(i);
%     Sol 8: CRLB
    S_all{i}(8,:) = alg_crlb_greedy(snr, D, M);
    
    % Sol 6: threshold greedy. single point
    S_all{i}(6,:) = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);
    
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S1', snr);
    MSE_avg_cell{i}(1) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
    MSE_avg_cell{i}(2) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
    MSE_avg_cell{i}(3) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S4', snr);
    MSE_avg_cell{i}(4) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S_all{i}(6,:)', snr);
    MSE_avg_cell{i}(6) = Es(2);
%     
    crlb = calculate_crlb(S_all{i}(8,:)', snr);
    MSE_avg_cell{i}(8) = crlb;
    
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution
% 
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim_cell{i}(1) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell{i}(2) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim_cell{i}(3) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim_cell{i}(4) = mean((results).^2);
%     
%     
    MSE_sim_cell{i}(8) = crlb;
    
    % overwrite...random antenna test..
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S_all{i}(2,:) = D(sort(randperm(21,M)));
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(2,:), snr);
    end
    MSE_sim_cell{i}(2) = mean((results).^2);
    
    
    % TRA-greedy-anchors = anchor_deg
%     results = zeros(1, sim_times);
%     for j = 1:sim_times
%         hat_u = sin(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
%         S_all{i}(5,:) = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr, D, M);
%         [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(5,:), snr);
%     end
%     MSE_sim_cell{i}(5) = mean((results).^2);
%     
% 
%   % TRA-greedy = anchor_deg
%     results = zeros(1, sim_times);
%     for j = 1:sim_times
%         hat_u = sin(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
%         S_all{i}(6,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
%         [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(6,:), snr);
%     end
%     MSE_sim_cell{i}(6) = mean((results).^2);
    
    
    
    % TRA-greedy = single anchor
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sin(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u*0;
        S_all{i}(7,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(7,:), snr);
    end
    MSE_sim_cell{i}(7) = mean((results).^2);

%     toc
end
toc
% save sim_M8.mat 
% load sim_M8_0804.mat 
% load sim_M8_1124.mat 

% load sim_los_M8_5K_1124.mat
% Approximation
% linewidth = 1.2;
% markersize = 8;
% MSE_avg = reshape(cell2mat(MSE_sim_cell), [9 length(snrdB)]);
% figure;semilogy(snrdB, MSE_avg(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_avg(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
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
% MSE_sim_cell = test = MSE_sim_cell;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [9 length(snrdB)]);

linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_sim(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '--oc', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('PSLC [\beta = 1]', 'PSLC [\beta = 0.75]', 'PSLC [\beta = 0.5]', ...
'ULA', 'TRA-G [e=5, n=3]', 'TRA-G2 [e=5, n=1]', 'TRA-G3 [e=1, n=1]', 'Best Case CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)

%% comparison
M = 6;
snrdB = -5:2.5:25; % M6

% Benchmark curves: 1-CRLB-PSL-1; 2-0.75; 3-0.5; 4-ULA; 5-TRA-Greedy; 
% 6-TRA-Greedy-p; 7-TRA-NN-p; 8-Optimal CRLB

MSE_avg_cell = zeros(8, length(snrdB));
MSE_sim_cell = zeros(8, length(snrdB));
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
    MSE_avg_cell(1,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
    MSE_avg_cell(2,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
    MSE_avg_cell(3,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S4', snr);
    MSE_avg_cell(4,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S6', snr);
    MSE_avg_cell(6,i) = Es(2);
    
    
    crlb = calculate_crlb(S8, snr);
    MSE_avg_cell(8,i) = crlb;
    
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim_cell(1,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell(2,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim_cell(3,i) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim_cell(4,i) = mean((results).^2);
    
    MSE_sim_cell(8,i) = crlb;


    % TRA-greedy-p
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S6 = alg_threshold_greedy(target_all(j), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S6, snr);
    end
%     [results, outputs] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell(6,i) = mean((results).^2);

    toc
end

% save simulation_M8.mat 
% load simulation_M8.mat 

% save simulation_M6.mat 
% load simulation_M6.mat 

% save sim_M6.mat 
% load sim_M6_0804.mat 
% load sim_M6_0915.mat
%% Approximation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_avg_cell(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_avg(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_avg_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_avg_cell(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('CRLB-PSL [\beta = 1]', 'CRLB-PSL [\beta = 0.75]', 'CRLB-PSL [\beta = 0.5]', ...
'ULA', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)
%% Simulation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_sim_cell(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_sim(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim_cell(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


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

MSE_avg_cell = zeros(8, length(snrdB));
MSE_sim_cell = zeros(8, length(snrdB));
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
    MSE_avg_cell(1,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
    MSE_avg_cell(2,i) = Es(2);

    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
    MSE_avg_cell(3,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S4', snr);
    MSE_avg_cell(4,i) = Es(2);
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S6', snr);
    MSE_avg_cell(6,i) = Es(2);
    
    crlb = calculate_crlb(S8, snr);
    MSE_avg_cell(8,i) = crlb;
    
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa)*180/pi); % uniform distribution

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim_cell(1,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell(2,i) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim_cell(3,i) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim_cell(4,i) = mean((results).^2);
    
    MSE_sim_cell(8,i) = crlb;


    % TRA-greedy-p
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S6 = alg_threshold_greedy(target_all(j), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S6, snr);
    end
%     [results, outputs] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell(6,i) = mean((results).^2);

    toc
end

% save simulation_M8.mat 
% load simulation_M8.mat 

% save simulation_M6.mat 
% load simulation_M6.mat 

% save sim_M4.mat 
% load sim_M6_0803.mat 
% load sim_M4_0915.mat
%% Approximation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_avg_cell(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_avg(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_avg_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_avg_cell(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('CRLB-PSL [\beta = 1]', 'CRLB-PSL [\beta = 0.75]', 'CRLB-PSL [\beta = 0.5]', ...
'ULA', 'TRA-Greedy', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)
%% Simulation
linewidth = 1.2;
markersize = 8;
figure;semilogy(snrdB, MSE_sim_cell(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
% hold on;semilogy(snrdB, MSE_sim(5,:), '-dr', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim_cell(7,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)

% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim_cell(6,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)


grid on;
legend('CRLB-PSL [\beta = 1]', 'CRLB-PSL [\beta = 0.75]', 'CRLB-PSL [\beta = 0.5]', ...
'ULA', 'TRA-Greedy-p', 'TRA-NN-p', 'CRLB')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)