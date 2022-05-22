close all;
clear all;
clc;
% random doa simulation
%% 
umax = 0.9; % TRA calculation range (-64 to 64 deg)
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
anchors = [-(tilde_u) -(tilde_u)/2 0 (tilde_u)/2 (tilde_u)]; % 5 anchors to add hat_u

thetam = floor(asind(umax));
sim_times = 5000;
rng(1);
t = datetime('now')

%% comparison

M = 6;
snrdB = [5 10 15];
snr1 = snrdB(1);
snr2 = snrdB(2);
snr3 = snrdB(3);

iter_vec = 1:8;

% Benchmark curves: 1-PSL-C-1.0; 2-0.85; 3-0.5; 4-ULA; 
% 5-TRA-Greedy [1 Anchor, perfect];  
% 6-TRA-Greedy [1 Anchor with error]; 
% 7-TRA-Greedy [5 Anchor with error];
% 8-Optimal CRLB; 
% 9-Luck; 
% 10-NN run with other part of code

Results_sim_cell = cell(1,sim_times);
target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % large angle area (fig-3)
% target_all = target_doa + (rand(1,sim_times)-0.5)*2*asind(0.1);     % a smaller angle candidate area


S_all = cell(1,length(iter_vec));
for i = 1:length(Results_sim_cell)
    Results_sim_cell{i} = zeros(10, length(iter_vec));
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


snr_error = 1;

tic
% for i = 1:sim_times
parfor i = 1:sim_times
    if mod(i, 10) == 0
        disp(i);
    end
    target = target_all(i);

    % ********* snr 1**************
    % Sol 1: 1-anchor inaccurate, tilde_u = 0.1;
    results = zeros(1, length(iter_vec));
    for j = 1:length(iter_vec)
        if j == 1
            S_cur = S4;
        else
            hat_u = mean(sind(results(1:j-1)));
%             hat_u = median(sind(results(1:j-1)));
            S_cur = alg_threshold_greedy(asind(hat_u), ref_doa, snr1, D, M);
        end
        [~, results(j), ~] = batch_doa_simulation(target, thetam, S_cur, snr1);
    end
    Results_sim_cell{i}(1,:) = results;
    
    % Sol 2: 5 anchors
    for j = 2:length(iter_vec)
        if j == 1
            S_cur = S4;
        else
            hat_u = mean(sind(results(1:j-1)));
%             hat_u = median(sind(results(1:j-1)));
            S_cur = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr1, D, M);
        end
        [~, results(j), ~] = batch_doa_simulation(target, thetam, S_cur, snr1);
    end
    Results_sim_cell{i}(2,:) = results;
    
     % ********* snr 2**************
    % Sol 1: 1-anchor inaccurate, tilde_u = 0.1;
    results = zeros(1, length(iter_vec));
    for j = 1:length(iter_vec)
        if j == 1
            S_cur = S4;
        else
            hat_u = mean(sind(results(1:j-1)));
%             hat_u = median(sind(results(1:j-1)));
            S_cur = alg_threshold_greedy(asind(hat_u), ref_doa, snr2, D, M);
        end
        [~, results(j), ~] = batch_doa_simulation(target, thetam, S_cur, snr2);
    end
    Results_sim_cell{i}(3,:) = results;
    
    % Sol 2: 5 anchors
    for j = 2:length(iter_vec)
        if j == 1
            S_cur = S4;
        else
            hat_u = mean(sind(results(1:j-1)));
%             hat_u = median(sind(results(1:j-1)));

            S_cur = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr2, D, M);
        end
        [~, results(j), ~] = batch_doa_simulation(target, thetam, S_cur, snr2);
    end
    Results_sim_cell{i}(4,:) = results;
    
    
     % ********* snr 3**************
    % Sol 1: 1-anchor inaccurate, tilde_u = 0.1;
    results = zeros(1, length(iter_vec));
    for j = 1:length(iter_vec)
        if j == 1
            S_cur = S4;
        else
            hat_u = mean(sind(results(1:j-1)));
%             hat_u = median(sind(results(1:j-1)));

            S_cur = alg_threshold_greedy(asind(hat_u), ref_doa, snr3, D, M);
        end
        [~, results(j), ~] = batch_doa_simulation(target, thetam, S_cur, snr3);
    end
    Results_sim_cell{i}(5,:) = results;
    
    % Sol 2: 5 anchors
    for j = 2:length(iter_vec)
        if j == 1
            S_cur = S4;
        else
            hat_u = mean(sind(results(1:j-1)));
%             hat_u = median(sind(results(1:j-1)));

            S_cur = alg_threshold_greedy(asind(hat_u + anchors), ref_doa, snr3, D, M);
        end
        [~, results(j), ~] = batch_doa_simulation(target, thetam, S_cur, snr3);
    end
    Results_sim_cell{i}(6,:) = results;
    
    
    
    % Sol 7: CRLB
    S_all{i}(7,:) = alg_crlb_greedy(snr1, D, M);
    % Best case CRLB
    Results_sim_cell{i}(7,:) = calculate_crlb(S_all{i}(7,:)', snr1);
        
    % Sol 8: CRLB
    S_all{i}(8,:) = alg_crlb_greedy(snr2, D, M);
    % Best case CRLB
    Results_sim_cell{i}(8,:) = calculate_crlb(S_all{i}(8,:)', snr2);
        
    % Sol 9: CRLB
    S_all{i}(9,:) = alg_crlb_greedy(snr3, D, M);
    % Best case CRLB
    Results_sim_cell{i}(9,:) = calculate_crlb(S_all{i}(9,:)', snr3);
        
    
end
toc
t = datetime('now')



%%
% save fig_7_500.mat
% save fig_7_5000.mat
% save fig_7_5000_small.mat

% load fig_7_5000.mat
% load fig_7_5000_small.mat

% Simulation
linewidth = 1.2;
markersize = 6;

Results_sim = zeros(sim_times, length(iter_vec), 10);
for i = 1:sim_times
    for j = 1:10
        Results_sim(i, :, j) = Results_sim_cell{i}(j,:);
    end
end

Errors_sim = zeros(10, length(iter_vec));
for i = 1:6
    Errors_sim(i,:) = mean(deg2rad(Results_sim(:,:,i) - target_all').^2);
end

Errors_sim(7,:) = Results_sim(1,:,7);
Errors_sim(8,:) = Results_sim(1,:,8);
Errors_sim(9,:) = Results_sim(1,:,9);


figure;semilogy(iter_vec, Errors_sim(1,:), '--rs', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(iter_vec, Errors_sim(2,:), '-rs', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(iter_vec, Errors_sim(7,:),'--k', 'Linewidth', linewidth, 'MarkerSize',markersize)

hold on;semilogy(iter_vec, Errors_sim(3,:), '--bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(iter_vec, Errors_sim(4,:), '-bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(iter_vec, Errors_sim(8,:),'-.k', 'Linewidth', linewidth, 'MarkerSize',markersize)

hold on;semilogy(iter_vec, Errors_sim(5,:), '--gd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(iter_vec, Errors_sim(6,:), '-gd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(iter_vec, Errors_sim(9,:),'-k', 'Linewidth', linewidth, 'MarkerSize',markersize)



legend('$N_A$=1 (5 dB)', '$N_A$=5 (5 dB)', 'CRLB (5 dB)', '$N_A$=1 (10 dB)', '$N_A$=5 (10 dB)', 'CRLB (10 dB)', ...
    '$N_A$=1 (15 dB)', '$N_A$=5 (15 dB)', 'CRLB (15 dB)', 'NumColumns',3, 'Interpreter','Latex');


grid on;
axis([1 8 10^-5.1 10^-1])
xlabel('Measurement Index', 'Interpreter','Latex');
ylabel('MSE', 'Interpreter','Latex');
grid on; set(gca,'FontSize',18);
set(gcf,'position',[100,100,600,350])

%%

