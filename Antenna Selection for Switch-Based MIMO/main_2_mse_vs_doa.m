close all;
clear all;
clc;
%%
global dv;
dv = [0 1 2 3 4 5]';   % in the unit of lambda/2

theta_max = 90;
snrdB   = -10:2.5:30; % different SNR for plotting

rng(1);
results_mat = zeros(7,length(snrdB)); % 1-pdp; 2-MLE; 3-MUSIC
time_all = zeros(1,size(results_mat,1));    % running time

% estimations
sim_times = 1000;
for i_snr = 1:length(snrdB)
    snr = snrdB(i_snr)
    % first DOA
    doa_center = 0;
    target_all = doa_center + (rand(1,sim_times)-0.5)*1; % uniform distribution
    results = batch_doa_simulation(target_all, theta_max, dv', snr);
    results_mat(1, i_snr) = mean((results).^2);
    
    doa_center = 30;
    target_all = doa_center + (rand(1,sim_times)-0.5)*1; % uniform distribution
    results = batch_doa_simulation(target_all, theta_max, dv', snr);
    results_mat(2, i_snr) = mean((results).^2);
    
    doa_center = 70;
    target_all = doa_center + (rand(1,sim_times)-0.5)*1; % uniform distribution
    results = batch_doa_simulation(target_all, theta_max, dv', snr);
    results_mat(3, i_snr) = mean((results).^2);
    
    % figure;plot(results)

end

%%

%% Estimated 
target_doa = 0;
ref_angle = -theta_max:1:theta_max;
ref_doa = sin(ref_angle/180*pi);
% 
% [E_old, ~, ~, ~, ~] = alg_threshold_region_old(target_doa, ref_doa, d, snrdB);
[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, dv, snrdB);

% figure;plot(ref_doa, V); hold on;
% plot(ref_doa(Vmax), V(Vmax),'ro')
% xlabel('sin(u)'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');

%%
% save data_fig_02.mat;
% load data_fig_02.mat;
%% CRLB
result = 0;
M = length(dv);
snr = M*10.^(snrdB/10);
for i = 1:M
    result = result + (dv(i)-mean(dv))^2;
end
result = result/M;

CRB = 1/2/pi^2/result./snr;
% crlb = CRB*(180/pi/cos(39.5/180*pi))^2;
crlb = CRB;
% hold on;semilogy(snrdB, crlb);
plot_ind = 1:length(CRBs);
linewidth = 1.5;
markersize = 8;

figure;
semilogy(snrdB(plot_ind),(results_mat(1, plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),(results_mat(2, plot_ind)),'-^r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),(results_mat(3, plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(results_mat(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(results_mat(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(rvec6(plot_ind)),'--^c','LineWidth',linewidth,'MarkerSize',markersize);
% semilogy(snrdB,(crlb),'-sk','LineWidth',linewidth,'MarkerSize',markersize);
% hold on; semilogy(snrdB(plot_ind), E_old(1, plot_ind),'--^c','LineWidth',2);
% hold on; semilogy(snrdB(plot_ind), E_old(2, plot_ind),'-dc','LineWidth',2);
% hold on; semilogy(snrdB(plot_ind), E_new(2, plot_ind),'--+r','LineWidth',2);
hold on; semilogy(snrdB(plot_ind), crlb(plot_ind),'--k','LineWidth',2);

% hold on;semilogy(snrdB(plot_ind),(results_mat(1, plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('MSE');
% legend('PDP','MLE (0.2^\circ)','MUSIC (1^\circ)', 'CLRB');
% legend('DOA = 0 [rad]', 'DOA = 0.5 [rad]', 'DOA = 1.2 [rad]', 'CLRB');
legend('u = 0 (\theta=0\circ)', 'u = 0.5 (\theta = 30\circ)', 'u = 0.94 (\theta=70\circ)', 'CLRB');

grid on; set(gca,'FontSize',12)
% set(gcf,'position',[100,100,420*1.2,400*1.2])
% set(gca,'FontSize',12)


% print -dpng -r600 fig-2-diffdoas.png
