% calculate full set, unique set and greedy set
% code for Table 1.
close all;
clear all;
clc;
%%
N = 11;
M = 4;
L = (N-1)/4;

Set_full = nchoosek(1:N, M);
K = size(Set_full,1);
index_all = zeros(K, N);

for i = 1:size(Set_full,1)
    index = zeros(1,N);
    index(Set_full(i,:)) = 1;
    index = alg_align_index(index);
    index_all(i,:) = index;
end

UN = unique(index_all,'row');
greedy = sum(M+1:N);

disp(num2str(N) + " choose " + num2str(M) + ":");
disp("Full set = " + num2str(K) + " | " + "Unique Set = " + num2str(size(UN,1)) ...
     + " | " + "Ratio = " + num2str(size(UN,1)/K, "%4.4f"));
disp("Greedy set = " + num2str(greedy));

UN_value = UN*[2.^(0:(N-1))]';
% out = find(UN_value == 131199)
% UN(out,:)
% UN(out,:)*[2.^(0:(N-1))]'

nchoosek(N, M);
unique_lut = UN;
value_lut = UN_value;

%% calculate the number of switches
tilde_M = zeros(M, N);
for i = 1:size(UN, 1)
    b_i = UN(i,:);
    ind = find(b_i==1);
    for j = 1:length(ind)
        tilde_M(j, ind(j)) = 1;
    end
end

disp("Number of switches: " + num2str(sum(tilde_M(:))) + " | " + "Full connections: " + num2str(M*N)); % number of switches
% disp(M*N);              % number of full connections
% sum(tilde_M(:))/(M*N);

%% algorithm to calculate unique sets.
aligned_set = alg_unique_calculation(N, M);
% nchoosek(N-1,M-1)
repeated_set = alg_redundant_calculation(N-2, M-2);
unique_set = aligned_set - repeated_set;
repeated_set = alg_redundant_calculation(N-2, M-2);
repeated_set;

rr = [];
for n = (M-1):(N-2)
    m = M-2;
    o_N = mod(n, 2);
    o_M = mod(m, 2);
    flag = (o_N == 0 && o_M == 1);
    rr = [rr (nchoosek(n, m) - (1-flag)*nchoosek(floor(n/2), floor(m/2)) )/2];
end
sum(rr);
unique_set = nchoosek(N-1, M-1) - sum(rr);
%% CRLB
% sidelobe_lut = zeros(1, size(UN,1));
% crlb_lut = zeros(size(UN,1), length(snrdB));
% target_doa = 0;
% ref_doa = -1:0.01:1;
% 
% sidelobe = zeros(size(UN));
% crlb = zeros(size(UN));
% 
% tic
% for i = 1:size(UN,1)
%     
%     array = UN(i,:);
%     S = D(array==1)';
%     
%     % calculate sidelobes
%     u0 = sin(target_doa/180*pi);    % target DOA in rad
%     a0 = exp(-1j*pi*S*u0);  % steering vector
%     tx = 1;         % transmitted signal, single snapshot
%     sigma = 0;
%     n = sigma*randn(size(a0)) + 1j*sigma*randn(size(a0));   % zero noise
%     x = a0*tx+n;
%     V = zeros(1,length(ref_doa));
%     for doa_i = 1:length(ref_doa)
%         u = ref_doa(doa_i);
%         au = exp(-1j*pi*S*u);
%         V(doa_i) = abs(au'*x)^2/(au'*au);
%     end
%     % figure;plot(V)
%     V2 = [0 V 0];   % pad zeros
%     diff1 = sign(diff(V2));
%     diff2 = diff(diff1);
%     Vmax = find(diff2==-2);
%     [B, I] = sort(V(Vmax),'descend');
%     if (length(B) == 1)
%         sidelobe_lut(i) = 0;
%     else
%         sidelobe_lut(i) = B(2)/B(1);
%     end
%     crlb_lut(i,:) = calculate_crlb(S, snrdB);
% end
% toc
% figure;plot(sidelobe_lut)
% 
% 
% % save lut_8_sensors_center.mat sidelobe_lut crlb_lut unique_lut value_lut snrdB
% % save lut_6_sensors_center.mat sidelobe_lut crlb_lut unique_lut value_lut snrdB
% % save lut_4_sensors_center.mat sidelobe_lut crlb_lut unique_lut value_lut snrdB
% 


