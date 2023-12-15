%% makeABbar_HW Function Test Code %%
clc
clear;

%% Initial Parameter %%
% Set the File Path
file_path =

% Load the .mat File
load(file_path);

tau = [3 6 9 12 18 24 30 36 60 120]';
k = 3;
GQ = [1	0	0; 0	0.935	0.0628; 0	0	0.935];

lambda = [-1;-1;-2];
Omega = [0.1664 0 0 ; 0	0.0272 0 ; 0 0	0.0064];
L = cholmod(Omega);

%% Function %%
[abar, bbar, Abar, Bbar] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, tau, k);
figure
plot(abar);
figure
plot(bbar);