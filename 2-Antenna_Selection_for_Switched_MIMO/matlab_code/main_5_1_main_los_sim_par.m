close all;
clear all;
clc;
% random doa simulation
%% 
umax = 0.9; % TRA calculation range
ref_doa = -umax:0.025:umax;
% ref_doa = -umax:0.01:umax;

L = 10;
interval = 0.5;
TRA_max = 0.9;
TRA_interval = 0.1;
D = 0:interval:L;

target_doa = 0;
max_doa = 0.9;  % simulation DOA max range, sind(65) = 0.9063
tilde_u = 0.1;   % DOA estimation error, 50~75 at 60 deg, -6~6 at 0 deg
small_u = tilde_u*0.2;   % compare the effect of DOA error
% anchors = [-(tilde_u) 0 (tilde_u)]*0.75; % 3 anchors to add hat_u
anchors = [-(tilde_u) -(tilde_u)/2 0 (tilde_u)/2 (tilde_u)]; % 3 anchors to add hat_u

thetam = floor(asind(umax));
sim_times = 5000;
rng(0);
t = datetime('now')

%% comparison
M = 4;
snrdB = 0:2.5:30;  % M4

% Benchmark curves: 1-PSL-C-1; 2-0.85; 3-0.5; 4-ULA; 
% 5-TRA-Greedy [N_A = 3]; 6-TRA-Greedy [N_A = 1]; % 7-TRA-Greedy [perfect]; 
% 8-Optimal CRLB; 9-Luck; 10-NN

MSE_sim_cell = cell(1,length(snrdB));
S_all = cell(1,length(snrdB));
for i = 1:length(MSE_sim_cell)
    MSE_sim_cell{i} = zeros(1, 10);
    S_all{i} = zeros(10, M);
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
% parfor i = 6:12
    snr = snrdB(i)
    
    % Sol 8: CRLB
    S_all{i}(8,:) = alg_crlb_greedy(snr, D, M);
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution

    % simulation of benchmarks (1-4)
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim_cell{i}(1) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell{i}(2) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim_cell{i}(3) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim_cell{i}(4) = mean((results).^2);
    

    % simulation of TRA-G (5-7)
    % 1-anchor perfect
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*small_u;
        S_all{i}(5,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(5,:), snr);
    end
    MSE_sim_cell{i}(5) = mean((results).^2);

    % 1-anchor inaccurate
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(6,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(6,:), snr);
    end
    MSE_sim_cell{i}(6) = mean((results).^2);

    % 3-anchor inaccurate
    results = zeros(1, sim_times);
    for j = 1:sim_times
%         hat_u = target_all(j) + (rand(1,1)-0.5)*2*5;
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(7,:) = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(7,:), snr);
    end
    MSE_sim_cell{i}(7) = mean((results).^2);
    
    % Best case CRLB
    MSE_sim_cell{i}(8) = calculate_crlb(S_all{i}(8,:)', snr);
        
    % luck...random antenna test..
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S_all{i}(9,:) = D(sort(randperm(21,M)));
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(9,:), snr);
    end
    MSE_sim_cell{i}(9) = mean((results).^2);
%     toc
end
toc
t = datetime('now')

% save sim_M4.mat

% load sim_M4.mat

% Simulation
linewidth = 1.2;
markersize = 6;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [10 length(snrdB)]);

figure;semilogy(snrdB, MSE_sim(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(5,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-.bs', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)

grid on;
legend('ULA', 'PSL-C [\delta = 1]', 'PSL-C [\delta = 0.85]', 'PSL-C [\delta = 0.5]', 'Random', ...
'Best Case CRLB', 'TRA-G [N=1 p]', 'TRA-G [N=1 e]', 'TRA-G [N=5 e]')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)

%% comparison
M = 6;
snrdB = -5:2.5:25;  % M8

% Benchmark curves: 1-PSL-C-1; 2-0.85; 3-0.5; 4-ULA; 
% 5-TRA-Greedy [N_A = 3]; 6-TRA-Greedy [N_A = 1]; % 7-TRA-Greedy [perfect]; 
% 8-Optimal CRLB; 9-Luck; 10-NN

MSE_sim_cell = cell(1,length(snrdB));
S_all = cell(1,length(snrdB));
for i = 1:length(MSE_sim_cell)
    MSE_sim_cell{i} = zeros(1, 10);
    S_all{i} = zeros(10, M);
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
% parfor i = 11

    snr = snrdB(i)
    
    % Sol 8: CRLB
    S_all{i}(8,:) = alg_crlb_greedy(snr, D, M);
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution

    % simulation of benchmarks (1-4)
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim_cell{i}(1) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell{i}(2) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim_cell{i}(3) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim_cell{i}(4) = mean((results).^2);
    

    % simulation of TRA-G (5-7)
    % 1-anchor perfect
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*small_u;
        S_all{i}(5,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(5,:), snr);
    end
    MSE_sim_cell{i}(5) = mean((results).^2);

    % 1-anchor inaccurate
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(6,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(6,:), snr);
    end
    MSE_sim_cell{i}(6) = mean((results).^2);

    % 3-anchor inaccurate
    results = zeros(1, sim_times);
    for j = 1:sim_times
%         hat_u = target_all(j) + (rand(1,1)-0.5)*2*5;
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(7,:) = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(7,:), snr);
    end
    MSE_sim_cell{i}(7) = mean((results).^2);
    
    % Best case CRLB
    MSE_sim_cell{i}(8) = calculate_crlb(S_all{i}(8,:)', snr);
        
    % luck...random antenna test..
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S_all{i}(9,:) = D(sort(randperm(21,M)));
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(9,:), snr);
    end
    MSE_sim_cell{i}(9) = mean((results).^2);
    
end
toc
t = datetime('now')


% save sim_M6.mat

% load sim_M6.mat

% Simulation
linewidth = 1.2;
markersize = 6;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [10 length(snrdB)]);

figure;semilogy(snrdB, MSE_sim(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(5,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-.bs', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)

grid on;
legend('ULA', 'PSL-C [\delta = 1]', 'PSL-C [\delta = 0.85]', 'PSL-C [\delta = 0.5]', 'Random', ...
'Best Case CRLB', 'TRA-G [N=1 p]', 'TRA-G [N=1 e]', 'TRA-G [N=5 e]')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)

%% comparison
M = 8;
snrdB = -5:2.5:25;  % M8

% Benchmark curves: 1-PSL-C-1; 2-0.85; 3-0.5; 4-ULA; 
% 5-TRA-Greedy [N_A = 3]; 6-TRA-Greedy [N_A = 1]; % 7-TRA-Greedy [perfect]; 
% 8-Optimal CRLB; 9-Luck; 10-NN

MSE_sim_cell = cell(1,length(snrdB));
S_all = cell(1,length(snrdB));
for i = 1:length(MSE_sim_cell)
    MSE_sim_cell{i} = zeros(1, 10);
    S_all{i} = zeros(10, M);
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
    snr = snrdB(i)
    
    % Sol 8: CRLB
    S_all{i}(8,:) = alg_crlb_greedy(snr, D, M);
    
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution

    % simulation of benchmarks (1-4)
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S1, snr);
    MSE_sim_cell{i}(1) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S2, snr);
    MSE_sim_cell{i}(2) = mean((results).^2);

    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S3, snr);
    MSE_sim_cell{i}(3) = mean((results).^2);
    
    [results, ~, ~] = batch_doa_simulation(target_all, thetam, S4, snr);
    MSE_sim_cell{i}(4) = mean((results).^2);
    

    % simulation of TRA-G (5-7)
    % 1-anchor perfect
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*small_u;
        S_all{i}(5,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(5,:), snr);
    end
    MSE_sim_cell{i}(5) = mean((results).^2);

    % 1-anchor inaccurate
    results = zeros(1, sim_times);
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(6,:) = alg_threshold_greedy(asind(hat_u), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(6,:), snr);
    end
    MSE_sim_cell{i}(6) = mean((results).^2);

    % 3-anchor inaccurate
    results = zeros(1, sim_times);
    for j = 1:sim_times
%         hat_u = target_all(j) + (rand(1,1)-0.5)*2*5;
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(7,:) = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(7,:), snr);
    end
    MSE_sim_cell{i}(7) = mean((results).^2);
    
    % Best case CRLB
    MSE_sim_cell{i}(8) = calculate_crlb(S_all{i}(8,:)', snr);
        
    % luck...random antenna test..
    results = zeros(1, sim_times);
    for j = 1:sim_times
        S_all{i}(9,:) = D(sort(randperm(21,M)));
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(9,:), snr);
    end
    MSE_sim_cell{i}(9) = mean((results).^2);
%     toc
end
toc
t = datetime('now')

% save sim_M8.mat 

% load sim_M8.mat

% Simulation
linewidth = 1.2;
markersize = 6;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [10 length(snrdB)]);

figure;semilogy(snrdB, MSE_sim(4,:), '-+g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--*g', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-or', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-->r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(5,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '--b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-.bs', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_avg(5,:),'->b', 'Linewidth', linewidth, 'MarkerSize',markersize)

grid on;
legend('ULA', 'PSL-C [\delta = 1]', 'PSL-C [\delta = 0.85]', 'PSL-C [\delta = 0.5]', 'Random', ...
'Best Case CRLB', 'TRA-G [N=1 p]', 'TRA-G [N=1 e]', 'TRA-G [N=5 e]')
xlabel('SNR (dB)')
ylabel('MSE')
grid on; set(gca,'FontSize',12)
