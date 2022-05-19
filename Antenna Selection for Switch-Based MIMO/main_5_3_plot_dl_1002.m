close all;
clear all;
clc;
% random doa simulation
%% 
% umax = 0.9;
% % ref_angle = -thetam:1:thetam;
% % ref_doa = sin(ref_angle/180*pi);
% ref_doa = -umax:0.01:umax;
% 
% thetam = floor(asind(umax));
% max_doa = 0.9;
% L = 10;
% interval = 0.5;
% TRA_max = 0.9;
% TRA_interval = 0.1;
% D = 0:interval:L;
% % snrdB = -40
% % crlb = calculate_crlb(S1, snr);
% 
% sim_times = 5000;
% rng(0);
%%
load sim_M4_1223.mat 
load model_continuous_M4_1123.mat; % NN model

tic
parfor i = 1:length(snrdB)
    snr = snrdB(i)
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
    results = zeros(1, sim_times);
    % NN
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(10,:) = alg_nn(abs(hat_u), (snr/30-0.5)*2, variables, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(10,:), snr);
    end
    MSE_sim_cell{i}(10) = mean((results).^2);
    
end
toc

% save sim_M4_1223.mat
%% Simulationlinewidth = 1.2;
% load sim_M4_1223.mat

linewidth = 1.5;
markersize = 10;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [10 length(snrdB)]);
figure;semilogy(snrdB, MSE_sim(4,:), '-*g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--sr', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-<r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-.+r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_sim(5,:), '--ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '-.b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(10,:),'-cd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k', 'Linewidth', linewidth, 'MarkerSize',markersize)

legends =get(gca,'Children');
axis([0 30 10^-8.4 1])
grid on; 
set(gca,'FontSize',18)
% legend([legend_children(1:10)], 'ULA', 'PSL-C [\delta = 1]', 'PSL-C [\delta = 0.85]', 'PSL-C [\delta = 0.5]', 'Random', ...
% 'TRA-G [N_A=1, \Deltau=0.02]', 'TRA-G [N_A=1, \Deltau=0.1]', 'TRA-G [N_A=5, \Deltau=0.1]', 'TRA-DL [\Deltau=0.1]',...
% 'Best Case CRLB', 'Location', 'SouthWest', 'NumColumns',2);
% xlabel('SNR (dB)')
% ylabel('MSE')
xlabel('SNR (dB)','Interpreter','Latex')
ylabel('MSE','Interpreter','Latex')

bounds = [27 28 10^-6.1 10^-5.7];
pos = [0.72 0.7 0.17 0.2];
vertex = [3 4];
p = gca;
% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
% y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1); 
% y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));

rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1], '-.k'); % Line to vertex 3
end
% legend hide
z = axes('position',pos);

box on % put box around new pair of axes
semilogy(snrdB, MSE_sim(4,:), '-*g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--sr', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-<r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-.>r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_sim(5,:), '--ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '-.b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(10,:),'-cd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

legend(legends(8:-1:1), 'ULA', 'PSL-C [$\delta$ = 1]', 'PSL-C [$\delta$ = 0.85]', 'PSL-C [$\delta$ = 0.5]', ...
'TRA-G [$N_A$=1]', 'TRA-G [$N_A$=5]', 'TRA-DL',...
'Best Case CRLB', 'Location', 'SouthWest', 'NumColumns',2, 'Interpreter','Latex');


% saveas(gcf,'sim_M4', 'epsc')
%%
load sim_M6_1223.mat 
load model_continuous_M6_1123.mat; % better

tic
parfor i = 1:length(snrdB)
    snr = snrdB(i)
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
    results = zeros(1, sim_times);
    % NN
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
%         hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*0;
        S_all{i}(10,:) = alg_nn(abs(hat_u), (snr/30-0.5)*2, variables, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(10,:), snr);
    end
    MSE_sim_cell{i}(10) = mean((results).^2);
    
end
toc
% save sim_M6_1223.mat

%% Simulationlinewidth = 1.2;
% load sim_M6_1223.mat

linewidth = 1.5;
markersize = 10;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [10 length(snrdB)]);
figure;semilogy(snrdB, MSE_sim(4,:), '-*g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--sr', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-<r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-.+r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_sim(5,:), '--ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '-.b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(10,:),'-cd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k', 'Linewidth', linewidth, 'MarkerSize',markersize)

grid on;
% legend('ULA', 'PSL-C [\delta = 1]', 'PSL-C [\delta = 0.85]', 'PSL-C [\delta = 0.5]', 'Random', ...
% 'TRA-G [N_A=1, \Deltau=0.02]', 'TRA-G [N_A=1, \Deltau=0.1]', 'TRA-G [N_A=5, \Deltau=0.1]', 'TRA-DL [\Deltau=0.1]',...
% 'Best Case CRLB', 'Location', 'SouthWest', 'NumColumns',2);
axis([-5 25 10^-8 1])
xlabel('SNR (dB)','Interpreter','Latex')
ylabel('MSE','Interpreter','Latex')

grid on; set(gca,'FontSize',18)
legends = get(gca,'Children');


bounds = [22 23 10^-5.7 10^-5.48];
pos = [0.72 0.7 0.17 0.2];
vertex = [3 4];
p = gca;
% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
% y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1); 
% y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));

rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1], '-.k'); % Line to vertex 3
end
% legend hide
z = axes('position',pos);

box on % put box around new pair of axes
semilogy(snrdB, MSE_sim(4,:), '-*g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--sr', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-<r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-.>r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_sim(5,:), '--ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '-.b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(10,:),'-cd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

legend(legends(8:-1:1), 'ULA', 'PSL-C [$\delta$ = 1]', 'PSL-C [$\delta$ = 0.85]', 'PSL-C [$\delta$ = 0.5]', ...
'TRA-G [$N_A$=1]', 'TRA-G [$N_A$=5]', 'TRA-DL',...
'Best Case CRLB', 'Location', 'SouthWest', 'NumColumns',2, 'Interpreter','Latex');

% print -dpng -r600 sim-M6.png


% saveas(gcf,'sim_M6', 'epsc')

%%
% load sim_M8.mat 
load sim_M8_1223.mat 
load model_continuous_M8_1123.mat;

tic
parfor i = 1:length(snrdB)
    snr = snrdB(i);
    % start simulation ************************************
    target_all = target_doa + (rand(1,sim_times)-0.5)*2*(asin(max_doa-tilde_u)*180/pi); % uniform distribution
    results = zeros(1, sim_times);
    % NN
    for j = 1:sim_times
        hat_u = sind(target_all(j)) + (rand(1,1)-0.5)*2*tilde_u;
        S_all{i}(10,:) = alg_nn(abs(hat_u), (snr/30-0.5)*2, variables, D, M);
        [results(j), output, ~] = batch_doa_simulation(target_all(j), thetam, S_all{i}(10,:), snr);
    end
    MSE_sim_cell{i}(10) = mean((results).^2);
end
toc

% save sim_M8_1223.mat

%% Simulationlinewidth = 1.2;
% load sim_M8_1223.mat

linewidth = 1.5;
markersize = 10;
MSE_sim = reshape(cell2mat(MSE_sim_cell), [10 length(snrdB)]);
figure;semilogy(snrdB, MSE_sim(4,:), '-*g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--sr', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-<r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-.+r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_sim(5,:), '--ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '-.bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-bs', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(10,:),'-cd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k', 'Linewidth', linewidth, 'MarkerSize',markersize)

grid on;
% legend('ULA', 'PSL-C [\delta = 1]', 'PSL-C [\delta = 0.85]', 'PSL-C [\delta = 0.5]', 'Random', ...
% 'TRA-G [N_A=1, \Deltau=0.02]', 'TRA-G [N_A=1, \Deltau=0.1]', 'TRA-G [N_A=5, \Deltau=0.1]', 'TRA-DL [\Deltau=0.1]',...
% 'Best Case CRLB', 'Location', 'SouthWest', 'NumColumns',2);
axis([-5 25 10^-8 1])
xlabel('SNR (dB)','Interpreter','Latex')
ylabel('MSE','Interpreter','Latex')

grid on; set(gca,'FontSize',18)

legends = get(gca,'Children');


bounds = [22 23 10^-5.8 10^-5.55];
pos = [0.72 0.7 0.17 0.2];
vertex = [3 4];
p = gca;
% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
% y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1); 
% y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));

rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1], '-.k'); % Line to vertex 3
end
% legend hide
z = axes('position',pos);

box on % put box around new pair of axes
semilogy(snrdB, MSE_sim(4,:), '-*g', 'Linewidth', linewidth, 'MarkerSize',markersize);
hold on;semilogy(snrdB, MSE_sim(1,:), '--sr', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(2,:), '-<r', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(3,:),'-.>r', 'Linewidth', linewidth, 'MarkerSize',markersize)
% hold on;semilogy(snrdB, MSE_sim(9,:), '--k', 'Linewidth', linewidth, 'MarkerSize',markersize+2);
% hold on;semilogy(snrdB, MSE_sim(5,:), '--ob', 'Linewidth', linewidth, 'MarkerSize',markersize-1)
hold on;semilogy(snrdB, MSE_sim(6,:), '-.b^', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(7,:), '-bo', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(10,:),'-cd', 'Linewidth', linewidth, 'MarkerSize',markersize)
hold on;semilogy(snrdB, MSE_sim(8,:),'-k+', 'Linewidth', linewidth, 'MarkerSize',markersize)
axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

legend(legends(8:-1:1), 'ULA', 'PSL-C [$\delta$ = 1]', 'PSL-C [$\delta$ = 0.85]', 'PSL-C [$\delta$ = 0.5]', ...
'TRA-G [$N_A$=1]', 'TRA-G [$N_A$=5]', 'TRA-DL',...
'Best Case CRLB', 'Location', 'SouthWest', 'NumColumns',2, 'Interpreter','Latex');

% print -dpng -r600 sim-M8.png

% saveas(gcf,'sim_M8', 'epsc')
