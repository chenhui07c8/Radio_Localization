close all% DOA estimation algorithms
clear all;
clc;
%% dv = [0 1 2 3 4 5]'; 
dv = [0 1 2 3 4 5]';   % in the unit of lambda/2
% S_lambda = [0 1 2 11 15 18];
% dv = [0 2 8 10]'
thetam = 90;

% snrdB   = 0:2.5:40; % different SNR for plotting
snrdB   = -15:2.5:25; % different SNR for plotting

results_mat = zeros(7,length(snrdB)); % 1-pdp; 2-MLE; 3-MUSIC
time_all = zeros(1,size(results_mat,1));    % running time

target_doa = 0;
ref_angle = -thetam:1:thetam;
ref_doa = sin(ref_angle/180*pi);

% [E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, dv, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}, {'u_2'}, {'u_3'}, {'u_1'}, {'u_4'}]
text(ref_doa(Vmax)-0.045, V(Vmax)+[-1 1 1 1 1]*0.6, xx, 'no info', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)

text(0.5, 5.3, '(K = 4)', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16);

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)

%% dv = [0 1 2 11 15 18];
dv = [0 1 2 11 15 18]';
target_doa = 0;
ref_angle = -thetam:1:thetam;
ref_doa = sin(ref_angle/180*pi);

% [E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, dv, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}]
text(ref_doa(Vmax(1))+0.05, V(Vmax(1))+[-1]*0.6, xx, 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)
text(0.5, 5.3, '(K = 14)', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16);

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)


%% for ploting
S_lambda = [0 1 2 3 4 5]*3;
target_doa = 0;
d = S_lambda';
ref_angle = -thetam:1:thetam;
ref_doa = sin(ref_angle/180*pi);

[E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, d, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}]
text(ref_doa(Vmax(1))+0.05, V(Vmax(1))+[-1]*0.6, xx, 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)
text(0.5, 5.3, '(K = 12)', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16);

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)

%%
S_lambda = [0 1 2 3 4 5]*2;
target_doa = 0;
d = S_lambda';
ref_angle = -thetam:1:thetam;
ref_doa = sin(ref_angle/180*pi);

[E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, d, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}]
text(ref_doa(Vmax(1))+0.05, V(Vmax(1))+[-1]*0.6, xx, 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)
text(0.5, 5.3, '(K = 8)', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16);

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)


%%
S_lambda = [0 3 6.5 9 12 15];
target_doa = 0;
d = S_lambda';
ref_angle = -thetam:1:thetam;
ref_doa = sin(ref_angle/180*pi);

[E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, d, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}]
text(ref_doa(Vmax(1))+0.05, V(Vmax(1))+[-1]*0.6, xx, 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)
text(0.5, 5.3, '(K = 14)', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16);

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)

%%
S_lambda = [0 2 6.5 9 12 15];
target_doa = 0;
d = S_lambda';
ref_angle = -thetam:1:thetam;
ref_doa = sin(ref_angle/180*pi);

[E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, d, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}]
text(ref_doa(Vmax(1))+0.05, V(Vmax(1))+[-1]*0.6, xx, 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)
text(0.5, 5.3, '(K = 12)', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16);

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)
set(gca,'Color','k')
set(gcf,'inverthardcopy','off'); 

