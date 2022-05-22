close all;
clear all;
clc;
% comparison between MSE simulation and approximation.

%%
global dv;
dv = [0 1 2 3 4 5]';   % antenna location vector [in lambda/2]
u_max = 1;   % max doa value, u = sin(theta)
theta_max = asin(u_max)*180/pi;  % max DOA value in degree
snrdB   = -15:2.5:25; % different SNR for plotting
errors_mat = zeros(1,length(snrdB)); % save DOA results using MLE

rng(0);
sim_times = 100;
for i_snr = 1:length(snrdB)
    snr = snrdB(i_snr)
    doa_center = 0; % DOA at 0
    target_all = doa_center + (rand(1,sim_times)-0.5)*2; % uniform distribution
    errors = batch_doa_simulation(target_all, theta_max, dv', snr);
    errors_mat(1, i_snr) = mean((errors).^2);
end
% figure;hist(target_all)
% figure;plot(errors_mat)
%%
% save data_fig_1.mat;
% load data_fig_1.mat;
%% Estimated 
target_doa = doa_center;
ref_angle = -theta_max:1:theta_max;
ref_doa = sin(ref_angle/180*pi);    % discrete sample on the DOA range to estimate expectation

[E_new, CRBs, threshold_snr, V, Vmax] = alg_threshold_region(target_doa, ref_doa, dv, snrdB);

figure;plot(ref_doa, V, 'Linewidth', 2); hold on;
plot(ref_doa(Vmax), V(Vmax),'ro','Linewidth', 2, 'MarkerSize', 10);
xx = [{'u_0'}, {'u_2'}, {'u_3'}, {'u_1'}, {'u_4'}]
text(ref_doa(Vmax)-0.045, V(Vmax)+[-1 1 1 1 1]*0.6, xx, 'no info', 'FontSize',12, ...
    'VerticalAlignment','middle', 'FontSize', 16)

xlabel('u'); ylabel('V(u)')
% figure;plot(ref_doa, V); hold on; plot(ref_doa(Vmax), V(Vmax),'ro')
% figure;semilogy(snrdB, Es(1,:),'-+c'); hold on; semilogy(snrdB, CRBs,'--k');
grid on; set(gca,'FontSize',16)
%%
% save fig_01.mat
% load fig_01.mat
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
plot_ind = 1:length(crlb);
linewidth = 1.5;
markersize = 6;

figure;
semilogy(snrdB(plot_ind),(errors_mat(1, plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(results_mat(2, plot_ind)),'-^r','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(results_mat(3, plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(results_mat(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(results_mat(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(rvec6(plot_ind)),'--^c','LineWidth',linewidth,'MarkerSize',markersize);
% semilogy(snrdB,(crlb),'-sk','LineWidth',linewidth,'MarkerSize',markersize);
hold on; semilogy(snrdB(plot_ind), E_new(2, plot_ind),'--^r','LineWidth',2);
% hold on; semilogy(snrdB(plot_ind), E_old(2, plot_ind),'-dc','LineWidth',2);
% hold on; semilogy(snrdB(plot_ind), E_new(2, plot_ind),'--+r','LineWidth',2);
hold on; semilogy(snrdB(plot_ind), crlb(plot_ind),'--k','LineWidth',2);

% hold on;semilogy(snrdB(plot_ind),(results_mat(1, plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);


fontsize = 12;
% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('MSE');
% legend('PDP','MLE (0.2^\circ)','MUSIC (1^\circ)', 'CLRB');
% legend('DOA = 0 [rad]', 'Approx-Original', 'Approx-Modified', 'CLRB');
axis([-15 25 10^-6 1]);
grid on; set(gca,'FontSize',12)
% set(gcf,'position',[100,100,420*1.2,400*1.2])
text(-12,10^-5, 'no info', 'VerticalAlignment','middle','FontSize',fontsize)
text(-2.2,10^-5, 'threshold', 'VerticalAlignment','middle','FontSize',fontsize)
text(10,10^-5, 'asymptotic', 'VerticalAlignment','middle','FontSize',fontsize)
% grid off;

hold on; semilogy([-2.5 -2.5], [1 10^-6],'-.k','LineWidth',0.5);
hold on; semilogy([5 5], [1 10^-6],'-.k','LineWidth',0.5);

legend('DOA MSE (u = 0)', 'Approximation', 'CLRB');
set(gca,'FontSize',12)

% print -dpng -r600 fig-1-regions.png