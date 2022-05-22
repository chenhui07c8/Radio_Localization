close all;
clear all;
clc;
% ULA sensor layout optimization through DOA threshold region performance approximation 
% ULA sensor layout optimization through DOA performance approximation and deep learning
% ULA DOA threshold region performance approximation and sensor layout optimization
%% 
angle_all = 45;
thetam = 90;
% ref_angle = -thetam:1:thetam;
% ref_doa = sin(ref_angle/180*pi);
snrdB = -5:2.5:20;
% threshold snr: 0:15
N = 21;
interval = 0.5;
L = (N-1)*interval;
um = 0.9;
ref_doa = -um:0.025:um;
D = 0:interval:L;
sim_times = 1000;
rng(0);
%% comparison
M = 4;
load model_continuous_M4_1123.mat

MSE_avg = zeros(5, length(snrdB));
target_doa = 0;
MSE_sim = zeros(5, length(snrdB));

TR_flag = 0;    % if set 1, run selection algorithm.
NN_flag = 0;

snrdB = 0; % M6
TR_limit = [-102.5 10020];    
NN_limit = [-102.5 10017.5];

times = zeros(5, sim_times);
snr = snrdB;

for i = 1:sim_times
    temp_D = randperm(21);
    d_temp = D(temp_D(1:M));
    tic;
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, d_temp', snrdB);
    times(1, i) = toc;
    
    tic;
    S3 = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);
    times(3, i) = toc;
    
    
    tic;
    S4 = alg_nn(sin(target_doa/180*pi), snr/30, variables, D, M);
    times(4, i) = toc;
end

% mc = nchoosek(N, M)
% mc = 7872
result = sum(times');

(N-M)*(N+M+1)/2;
% result(1)*mc*0.1*10
% result(2)*10
disp('***************  M=4  ***************');
% result(3)
% result(4)

% cvx
umax = 0.9;
M = 6;
time_cvx = zeros(1,sim_times);
delta = 0.5;
for i = 1:length(time_cvx)
    delta = 1;
%     [P_am2, x_am] = linear_amb(D,M,delta);
    tic;
    [S1, ~] = linear_iso_search(D,M,delta);
    time_cvx(i) = toc;
end
% sum(time_cvx)

result_M4 = [sum(time_cvx) result(3) result(4)]
%% comparison
M = 6;
load model_continuous_M6_1123.mat

MSE_avg = zeros(5, length(snrdB));
target_doa = 0;
MSE_sim = zeros(5, length(snrdB));

TR_flag = 0;    % if set 1, run selection algorithm.
NN_flag = 0;

snrdB = 0; % M6
TR_limit = [-102.5 10020];    
NN_limit = [-102.5 10017.5];

times = zeros(5, sim_times);
snr = snrdB;

for i = 1:sim_times
    temp_D = randperm(21);
    d_temp = D(temp_D(1:M));
    tic;
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, d_temp', snrdB);
    times(1, i) = toc;
    
    tic;
    S3 = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);
    times(3, i) = toc;
    
    
    tic;
    S4 = alg_nn(sin(target_doa/180*pi), snr/30, variables, D, M);
    times(4, i) = toc;
end

% mc = nchoosek(N, M)
% mc = 7872
result = sum(times');

(N-M)*(N+M+1)/2;
disp('***************  M=6  ***************');
% result(3)
% result(4)

% cvx
umax = 0.9;
M = 6;
time_cvx = zeros(1,sim_times);
delta = 0.5;
for i = 1:length(time_cvx)
    delta = 1;
%     [P_am2, x_am] = linear_amb(D,M,delta);
    tic;
    [S1, ~] = linear_iso_search(D,M,delta);
    time_cvx(i) = toc;
end
% sum(time_cvx)

result_M6 = [sum(time_cvx) result(3) result(4)]

%% comparison
M = 8;
load model_continuous_M8_1123.mat

MSE_avg = zeros(5, length(snrdB));
target_doa = 0;
MSE_sim = zeros(5, length(snrdB));

TR_flag = 0;    % if set 1, run selection algorithm.
NN_flag = 0;

snrdB = 0; % M6
TR_limit = [-102.5 10020];    
NN_limit = [-102.5 10017.5];

times = zeros(5, sim_times);
snr = snrdB;

for i = 1:sim_times
    temp_D = randperm(21);
    d_temp = D(temp_D(1:M));
    tic;
    
    [Es, ~, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, d_temp', snrdB);
    times(1, i) = toc;
    
    tic;
    S3 = alg_threshold_greedy(target_doa, ref_doa, snr, D, M);
    times(3, i) = toc;
    
    
    tic;
    S4 = alg_nn(sin(target_doa/180*pi), snr/30, variables, D, M);
    times(4, i) = toc;
end

% mc = nchoosek(N, M)
% mc = 7872
result = sum(times');

(N-M)*(N+M+1)/2;
% result(1)*mc*0.1*10
% result(2)*10
disp('***************  M=8  ***************');
% result(3)
% result(4)

% cvx
umax = 0.9;
M = 6;
time_cvx = zeros(1,sim_times);
delta = 0.5;
for i = 1:length(time_cvx)
    delta = 1;
%     [P_am2, x_am] = linear_amb(D,M,delta);
    tic;
    [S1, ~] = linear_iso_search(D,M,delta);
    time_cvx(i) = toc;
end
% sum(time_cvx)
result_M8 = [sum(time_cvx) result(3) result(4)]

% save data_fig_11.mat
% load data_fig_11.mat
