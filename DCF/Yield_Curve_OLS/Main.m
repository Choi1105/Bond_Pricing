%% Nelson-Siegel Model

% Y = Lambda * Beta + e

clear; 
clc; 

%% Step 1 : Make Lambda Matrix

lambda = 0.0609;
tau = [3 6 9 12 18 24 30 36 60 120];
[Lambda_M] = makeLambda(lambda,tau);

%% Step 2 : Estimate

Data = readmatrix("KOR_YC_2021_08_M.xlsx", "Range", "B2:K249");
Beta_Hat = inv(Lambda_M'*Lambda_M)*Lambda_M'*Data';
Beta = Beta_Hat';


%% PCA Part

%% Step 1: Principal Component Analysis
ym = readmatrix("KOR_YC_2021_08_M.xlsx", "Range","b2:k249");

% Data
Y = ym;
k = cols(Y);

% Standardization
stdY = standdc(Y);

% Sample correlation matrix
corrm = corrcoef(stdY);

% eigenvalue(D) and eigen vector(V)
[V,D] = eig(corrm);
  
% sorting
eigen_val = diag(D); % k by 1
[eigen_val, index] = sort(eigen_val, 'descend'); 
Vm = V(:, index); 

% compute PCA
PCm = stdY*Vm; % T by k

% Proportion 
cum_eigv = cumsum(eigen_val)/k;
disp('======================'); 
disp('Proportion of PCAs');
disp('----------------------');
disp(cum_eigv);
disp('======================');

%% Step 2: Plot
% Level, Slope, Curvature
Level = Y(:,end);
Slope = Y(:,1) - Y(:,end);
Curvature = (2*Y(:,5) - Y(:,end) - Y(:,1))/2;

figure
subplot(3, 1, 1)
plot([Level PCm(:,1) Beta(:,1)]);
legend('level', 'PCA', 'OLS');
title('Level with Estimation Plot')

subplot(3, 1, 2)
plot([Slope PCm(:,2) Beta(:,2)]);
legend('Slope', 'PCA', 'OLS');
title('Slope with Estimation Plot')

subplot(3, 1, 3)
plot([Curvature PCm(:,3) Beta(:,3)]);
legend('Curvature', 'PCA', 'OLS');
title('Curvature with Estimation Plot')

