

close all;
clear all;
clc;
%%
n1 = 2;
n2 = 3;
nm = 4;

Jm = randn(n1+n2,nm);

A = randn(n1)+2;
B = randn(n1, n2)+2;
C = transpose(B);
D = randn(n2);

J = [A B; C D];

inv1 = J^(-1);
inv1(1,1)


EJ = A-B*D^(-1)*C;
inv1 = J^(-1);
inv1(1,1)


%%


