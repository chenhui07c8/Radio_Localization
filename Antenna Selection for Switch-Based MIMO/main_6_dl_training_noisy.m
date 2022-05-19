% Create Training Dataset for TensorFlow
% The trained model can be found in model_continuous_M4/M6/M8
% alg_nn will run the neural network in Matlab
close all;
clear all;
clc;
%% Parameters
ref_u = -0.9:0.025:0.9;
ref_doa = ref_u;
snrdB = -5:2.5:40;
% threshold snr: 0:15
L = 10;
M = 4;
interval = 0.5;
D = 0:interval:L;
max_doa = 0.9;
% snrdB = -40
sim_times = 11000;
snr_min = 0;
snr_max = 30;
tilde_u = 0.1;   % DOA estimation error, 50~75 at 60 deg, -6~6 at 0 deg
anchors = [-(tilde_u) -(tilde_u)/2 0 (tilde_u)/2 (tilde_u)]; % 3 anchors to add hat_u


%% Generate training data: M4 discrete, 5-25 dB
doa_all = (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
M = 4;
% discrete training data
snr_all = randi([snr_min,snr_max],1,sim_times)*1; % 0 to 20 dB

% continuous training data.
% snr_all = rand(1,sim_times)*20; % 0 to 20 dB
% snr_all = rand(1,sim_times)*40-10; % 0 to 20 dB

index_all = zeros(sim_times, length(D));

tic
for i = 1:length(doa_all)
    if(mod(i,100) == 0)
        disp(i);
    end
    doa = doa_all(i);
%     doa = 0 + asin(-max_doa:0.2:max_doa)*180/pi;
    snr = snr_all(i);
    [S2, IND2, ~] = alg_threshold_greedy(asind(doa + anchors), ref_doa, snr, D, M);
    ind = alg_align_index(IND2);
    index_all(i,:) = ind;
end
toc
save data_N21_M4_discrete_10K.mat doa_all snr_all index_all

% figure;bar(sum(index_all)/1000);
% ylabel('Probability of Selection'); xlabel('Antenna Index')
% grid on; set(gca,'FontSize',12);
% A = unique(index_all, 'rows');
% size(A)
%% Generate training data: M4 continuous
doa_all = (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
M = 4;
% discrete training data
% snr_all = randi([0,20],1,sim_times)*1; % 0 to 20 dB

% continuous training data.
snr_all = snr_min + rand(1,sim_times)*(snr_max-snr_min); % 0 to 20 dB
% snr_all = rand(1,sim_times)*40-10; % 0 to 20 dB

index_all = zeros(sim_times, length(D));

tic
for i = 1:length(doa_all)
    if(mod(i,100) == 0)
        disp(i);
    end
    doa = doa_all(i);
%     doa = 0 + asin(-max_doa:0.2:max_doa)*180/pi;
    snr = snr_all(i);
    [S2, IND2, ~] = alg_threshold_greedy(asind(doa+anchors), ref_doa, snr, D, M);
    ind = alg_align_index(IND2);
    index_all(i,:) = ind;
end
toc
save data_N21_M4_continuous_10K.mat doa_all snr_all index_all



%% Generate training data: M6 discrete, 
doa_all = (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
M = 6;
% discrete training data
snr_all = randi([snr_min, snr_max],1,sim_times)*1; % 0 to 20 dB

% continuous training data.
% snr_all = rand(1,sim_times)*20; % 0 to 20 dB
% snr_all = rand(1,sim_times)*40-10; % 0 to 20 dB

index_all = zeros(sim_times, length(D));

tic
for i = 1:length(doa_all)
    if(mod(i,100) == 0)
        disp(i);
    end
    doa = doa_all(i);
%     doa = 0 + asin(-max_doa:0.2:max_doa)*180/pi;
    snr = snr_all(i);
    [S2, IND2, ~] = alg_threshold_greedy(asind(doa + anchors), ref_doa, snr, D, M);
    ind = alg_align_index(IND2);
    index_all(i,:) = ind;
end
toc
save data_N21_M6_discrete_10K.mat doa_all snr_all index_all

%% Generate training data: M6 discrete, 2-22 dB
doa_all = (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
M = 6;
% discrete training data
% snr_all = randi([2,22],1,sim_times)*1; % 0 to 20 dB

% continuous training data.
snr_all = snr_min + rand(1,sim_times)*(snr_max-snr_min); % 0 to 20 dB
% snr_all = rand(1,sim_times)*40-10; % 0 to 20 dB

index_all = zeros(sim_times, length(D));

tic
for i = 1:length(doa_all)
    if(mod(i,100) == 0)
        disp(i);
    end
    doa = doa_all(i);
%     doa = 0 + asin(-max_doa:0.2:max_doa)*180/pi;
    snr = snr_all(i);
    [S2, IND2, ~] = alg_threshold_greedy(asind(doa + anchors), ref_doa, snr, D, M);
    ind = alg_align_index(IND2);
    index_all(i,:) = ind;
end
toc
save data_N21_M6_continuous_10K.mat doa_all snr_all index_all

%%





%% Generate training data: M8 discrete, 
doa_all = (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
M = 8;
% discrete training data
snr_all = randi([snr_min, snr_max],1,sim_times)*1; % 0 to 20 dB

% continuous training data.
% snr_all = rand(1,sim_times)*20; % 0 to 20 dB
% snr_all = rand(1,sim_times)*40-10; % 0 to 20 dB

index_all = zeros(sim_times, length(D));

tic
for i = 1:length(doa_all)
    if(mod(i,100) == 0)
        disp(i);
    end
    doa = doa_all(i);
%     doa = 0 + asin(-max_doa:0.2:max_doa)*180/pi;
    snr = snr_all(i);
    [S2, IND2, ~] = alg_threshold_greedy(asind(doa + anchors), ref_doa, snr, D, M);
    ind = alg_align_index(IND2);
    index_all(i,:) = ind;
end
toc
save data_N21_M8_discrete_10K.mat doa_all snr_all index_all

%% Generate training data: M6 discrete, 2-22 dB
doa_all = (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
M = 8;
% discrete training data
% snr_all = randi([2,22],1,sim_times)*1; % 0 to 20 dB

% continuous training data.
snr_all = snr_min + rand(1,sim_times)*(snr_max-snr_min); % 0 to 20 dB
% snr_all = rand(1,sim_times)*40-10; % 0 to 20 dB

index_all = zeros(sim_times, length(D));

tic
for i = 1:length(doa_all)
    if(mod(i,100) == 0)
        disp(i);
    end
    doa = doa_all(i);
%     doa = 0 + asin(-max_doa:0.2:max_doa)*180/pi;
    snr = snr_all(i);
    [S2, IND2, ~] = alg_threshold_greedy(asind(doa + anchors), ref_doa, snr, D, M);
    ind = alg_align_index(IND2);
    index_all(i,:) = ind;
end
toc
save data_N21_M8_continuous_10K.mat doa_all snr_all index_all


%%