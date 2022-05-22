close all;
clear all;
clc;

%% 
angle_all = 0;      % DOA angle
um = 0.99;
ref_u = -um:0.01:um;
ref_doa = ref_u;
snrdB = -10:5:30;

L = 10;
N = 6;
interval = 0.5;
D = 0:interval:L;
%% quick comparison
MSE_avg = zeros(3, length(snrdB));
target_doa = angle_all(1);

for i = 1:length(snrdB)
    snr = snrdB(i);
    % sol 1: benchmark
%     S1 = D(1:N);
    S1 = 1:N;
    % sol 2: threshold region greedy
    S2 = alg_threshold_greedy(target_doa, ref_doa, snr, D, N);
    
    S3 = alg_crlb_greedy(snr, D, N);
    % sol 3: CRLB
    
    [Es, CRBs, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S1', snr);
%     Es(2)
    MSE_avg(1,i) = min(1, Es(2));
%     MSE_avg(3,i) = CRBs;

    [Es, CRBs, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S2', snr);
%     Es(2)
    MSE_avg(2,i) = min(1, Es(2));

    [Es, CRBs, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S3', snr);
%     Es(2)
    MSE_avg(3,i) = min(1, Es(2));


    crlb = calculate_crlb(S3, snr);
    MSE_avg(4,i) = crlb;
    
    crlb_ula = calculate_crlb(S1, snr);
    MSE_avg(5,i) = crlb_ula;
end

linewidth = 1.2;
markersize = 6;
figure;semilogy(snrdB, MSE_avg(1,:), '-ob', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_avg(2,:), '-^r', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_avg(3,:), '-dg', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(4,:),'-k', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_avg(5,:),'--k', 'Linewidth', linewidth, 'MarkerSize',markersize)

grid on;
legend('ULA', 'TR Greedy', 'CRLB Greedy', 'CRLB', 'CRLB-ULA')
xlabel('SNR [dB]');
ylabel('MSE');
grid on; set(gca,'FontSize',12)
%% optimal and greedy
L = 10;
N = 6;
interval = 0.5;
D = 0:interval:L;
snrdB = 30;
target_doa = 0;
ref_doa = -1:0.01:1;
IOU = zeros(1,5);   % greedy, 0.8 mc, 0.5 mc
Approx = zeros(1,5);

[S_g, Ind_g, ~] = alg_threshold_greedy(target_doa, ref_doa, snr, D, N);
[Es, CRBs, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, S_g', snr);
Val_g = Es(2);

Set_full = nchoosek(1:length(D), N);
Val_full = zeros(1, size(Set_full,1));
for i = 1:length(Val_full)
    [Es, CRBs, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, D(Set_full(i,:))', snr);
    Val_full(i) = Es(2);
end

[B, I] = sort(Val_full);
% figure;plot(Val_full);
min(Val_full);
% Val_g
temp = zeros(1,length(D));
temp(Set_full(I(1),:)) = 1;
Val_opt = B(1);
Ind_opt = temp;
% temp
% Ind_g
IOU(1) = sum(temp&Ind_g)/sum(temp|Ind_g);
Approx(1) = Val_g - min(Val_full);

% Monte Carlo 0.8
sim_times = 1000;
ratio = 0.75
select_times = round(ratio*length(Val_full));
Val_mc = zeros(1,sim_times);
IOU_mc = zeros(1,sim_times);
for sim_i = 1:sim_times 
    p = randperm(length(Val_full));
    [B, I] = sort(Val_full(p(1:select_times)));
    temp = zeros(1,length(D));
    temp(Set_full(p(I(1)),:)) = 1;
    IOU_mc(sim_i) = sum(temp&Ind_opt)/sum(temp|Ind_opt);
    Val_mc(sim_i) = B(1) - Val_opt;
end
IOU(2) = mean(IOU_mc);
Approx(2) = mean(Val_mc);

% MC 0.5
ratio = 0.5
select_times = round(ratio*length(Val_full));
Val_mc = zeros(1,sim_times);
IOU_mc = zeros(1,sim_times);
for sim_i = 1:sim_times 
    p = randperm(length(Val_full));
    [B, I] = sort(Val_full(p(1:select_times)));
    temp = zeros(1,length(D));
    temp(Set_full(p(I(1)),:)) = 1;
    IOU_mc(sim_i) = sum(temp&Ind_opt)/sum(temp|Ind_opt);
    Val_mc(sim_i) = B(1) - Val_opt;
end
IOU(3) = mean(IOU_mc);
Approx(3) = mean(Val_mc);

% MC 0.25
ratio = 0.25
select_times = round(ratio*length(Val_full));
Val_mc = zeros(1,sim_times);
IOU_mc = zeros(1,sim_times);
for sim_i = 1:sim_times 
    p = randperm(length(Val_full));
    [B, I] = sort(Val_full(p(1:select_times)));
    temp = zeros(1,length(D));
    temp(Set_full(p(I(1)),:)) = 1;
    IOU_mc(sim_i) = sum(temp&Ind_opt)/sum(temp|Ind_opt);
    Val_mc(sim_i) = B(1) - Val_opt;
end
IOU(4) = mean(IOU_mc);
Approx(4) = mean(Val_mc);


IOU
Approx+Val_opt

%% visualize the results
L = 10;
N = 6;
interval = 0.5;
D = 0:interval:L;
snr = 0;
% snr = 30;
IOU = zeros(1,5);   % greedy, 0.8 mc, 0.5 mc
Approx = zeros(1,5);
target_doa = 0;
ref_doa = -1:0.01:1;
Set_full = nchoosek(1:length(D), N);
Val_full = zeros(1, size(Set_full,1));
Aligned_full = zeros(length(Val_full), length(D));

for i = 1:length(Val_full)
    [Es, CRBs, ~, ~, ~] = alg_threshold_region(target_doa, ref_doa, D(Set_full(i,:))', snr);
    Val_full(i) = Es(2);
    temp = zeros(1, length(D));
    temp(Set_full(i, :)) = 1;
    Aligned_full(i, :) = alg_align_index(temp);
end

% temp = [1 1 0 1 0 0 0 0 1 1 1]
    
[B, I] = sort(Val_full);
min(Val_full)
B(1:5)

Aligned_full(I(1:5),:)


j = 1;
prev = zeros(1,length(D));
figure;
for i = 1:10
    now = Aligned_full(I(j),:);
    j=j+1;
    while(sum(prev&now)==N)
        j = j+1;
        now = Aligned_full(I(j),:);
    end
    prev = now;
    disp(now);
    plot([0 L], [i, i],'--b','Linewidth', 1.5, 'Markersize', 10);
    hold on;
    plot(D(logical(now)), i*ones(1, N),'rx','Linewidth', 1.5, 'Markersize', 10);
    hold on;
end

xlabel('Sensor layout [\lambda/2]')
ylabel('Layout index')
grid off; set(gca,'FontSize',12)

% figure;plot(C(1:10),'-o', 'Linewidth',2)
% xlabel('Top-10 layouts')
% ylabel('Approximated MSE');

%%