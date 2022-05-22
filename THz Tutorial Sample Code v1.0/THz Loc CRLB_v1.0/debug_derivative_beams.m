close all
clear all
clc;

%%
theta = -90:1:90;
pattern_mat = zeros(5, length(theta));
delta = 1;

% a = exp(-1j*2*pi*d/lambda)
theta0 = 0;
S = 0:1:5;    % in [lambda/2]
S = delta*S(1:6);
a = exp(1j*pi*(S')*sind(theta))/sqrt(length(S));
b = exp(1j*pi*S'*sind(theta0));
% figure;plot(theta, abs(a*b))
pattern_mat(1, :) = abs(a'*b);

S = 0:1:5;
S = delta*S(1:6);
a = exp(1j*pi*(S')*cosd(theta))/sqrt(length(S));
b = exp(1j*pi*S'*sind(theta0));
% hold on;plot(theta, abs(a*b))
pattern_mat(2, :) = abs(a'*b);

figure;semilogy(pattern_mat')

% S = 0:0.5:5;
% S = delta*S(1:11);
% a = exp(1j*pi*(S')*sind(theta))'/sqrt(length(S));
% b = exp(1j*pi*S'*sind(theta0));
% pattern_mat(3, :) = abs(a*b);

% S = [1:5:25 51:5:75];
% S = 0:5:25;
% S = delta*S(1:6);
% a = exp(1j*pi*(S')*sind(theta))'/sqrt(length(S));
% b = exp(1j*pi*S'*sind(theta0));
% pattern_mat(4, :) = abs(a*b);

% S = [1:5:25 51:5:75];
% S = 0:1:25;
% delta = 0.25;
% S = delta*S;
% a = exp(1j*pi*(S')*sind(theta))'/sqrt(length(S));
% % b = exp(1j*pi*S'*sind(theta0));
% b1 = exp(1j*pi*S(1:13)'*sind(30));
% b2 = exp(1j*pi*S(14:end)'*sind(-45));
% b = [b1; b2];
% pattern_mat(5, :) = abs(a*b);


%%
% colors = ["#003ADE", "#FF1F5B", "#00CD6C", "#AF58BA", "#F28522", "#FFC61E", ];
colors = ["#003ADE", "#FF1F5B", "#00CD6C", "#7F00FF", "#FFBF00", "#FFC61E", ];

% figure;
% polarplot(deg2rad(theta),pattern_mat(1,:), '-', 'color', colors(1), 'Linewidth', 1.4);
% hold on;polarplot(deg2rad(theta),pattern_mat(2,:), '--', 'color', colors(2), 'Linewidth', 1.1);
% hold on;polarplot(deg2rad(theta),pattern_mat(3,:), '-.', 'color', colors(3), 'Linewidth', 1.1);
% hold on;polarplot(deg2rad(theta),pattern_mat(4,:), '-', 'color', colors(4), 'Linewidth', 1.1);

figure;
polarplot(deg2rad(theta),pattern_mat(1,:), 'b-', 'Linewidth', 1.4);
hold on;polarplot(deg2rad(theta),pattern_mat(2,:), 'r--', 'Linewidth', 1.1);
hold on;polarplot(deg2rad(theta),pattern_mat(3,:), 'g-.', 'Linewidth', 1.1);
hold on;polarplot(deg2rad(theta),pattern_mat(4,:), '-', 'color', colors(5), 'Linewidth', 1.1);


% figure;plot(theta, pattern_mat)
thetalim([-90 90])

pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';

legend( '{$\mathring\Delta$=$0.5{\lambda}$/$\mathring{N}_A$=6/$D$=2.5$\lambda$}',...
        '{$\mathring\Delta$=$0.5{\lambda}$/$\mathring{N}_A$=2/$D$=0.5$\lambda$}', ...
        '{$\mathring\Delta$=$0.25{\lambda}$/$\mathring{N}_A$=11/$D$=2.5$\lambda$}', ...
        '{$\mathring\Delta$=$ 2.5{\lambda}$\hspace{1.3mm}/$\mathring{N}_A$= 6\hspace{1.3mm}/$D$=12.5$\lambda$}', ...
        'Location', 'South', 'Interpreter','Latex','NumColumns',2);

% Ambiguity can be solved by using a directional antenna and the sidelobes
% will be removed.

% xlabel('Number of RIS Elements');
% ylabel('OEB [Deg]');
% axis([-90 90 10^-0.6 10^1.6])

% grid on;
% set(gca,'fontsize', 14);
% set(gcf,'position', [100,100, 400*1.69, 400])

pax.FontSize = 14;
% set(gcf,'position',[100,100,600,600])
set(gcf,'position',[100,100,600,400])

% print -dpng -r600 sim_5_beamwidth.png

%%
N1 = 1;
N2 = 5;

theta1 = 2*asind(2.782/pi/N1)



N1 = 4
theta1 = 2*asind(2.782/pi/N1)

N2 = 20
theta2 = 2*asind(2.782/pi/N2)
%%
theta = -90:1:90;
pattern_mat = zeros(5, length(theta));
delta = 1;

% a = exp(-1j*2*pi*d/lambda)
S = 0:1:5;    % in [lambda/2]
S = delta*S(1:6);
a = exp(1j*pi*(S')*sind(theta))'/sqrt(length(S));

theta0 = -20;
b0 = exp(1j*pi*S'*sind(theta0));
theta1 = 40;
b1 = exp(1j*pi*S'*sind(theta1));
b = b0 + b1;
b = b./abs(b);

b = exp(1j*2*pi*rand(6,1));
figure;plot(theta, abs(a*b))

figure;plot(theta, abs(a*b0 + a*b1))


figure;plot(theta, abs(a*b + a*b) )

figure;plot(theta, abs(a*b0 + a*b1))


% figure;plot(theta, abs(a*b0))
% figure;plot(theta, abs(a*b1))

% pattern_mat(1, :) = abs(a*b);
%%
theta = -90:1:90;
pattern_mat = zeros(5, length(theta));
delta = 1;
N = 10;
S = 0:1:(N-1);    % in [lambda/2]
S = delta*S(1:N);
a = exp(1j*pi*(S')*sind(theta))'/sqrt(length(S));
b = exp(1j*2*pi*rand(N,1));
fftMat = fft(eye(N));
b = fftMat(:,4);
b = b./abs(b);
figure;plot(theta, abs(a*b))

%%
theta = -90:1:90;
pattern_mat = zeros(5, length(theta));
delta = 0.2;

% a = exp(-1j*2*pi*d/lambda)
theta0 = 0;
S = 0:5:25;    % in [lambda/2]
S = delta*S(1:6);
a = exp(1j*pi*(S')*sind(theta))'/sqrt(length(S));
b = exp(1j*pi*S'*sind(theta0));
% figure;plot(theta, abs(a*b))
pattern_mat(1, :) = abs(a*b);
