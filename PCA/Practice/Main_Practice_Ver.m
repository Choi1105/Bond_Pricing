

%% Principal Component Analysis
% X(t) = mu + phi*X(t-1) + v(t)
% where v(t) ~ N(0,sig2v)

% Y1(t) = gam1*X(t) + e1(t)
% Y2(t) = gam2*X(t) + e2(t)
% Y3(t) = gam3*X(t) + e3(t)
% Y4(t) = gam4*X(t) + e4(t)
% Y5(t) = gam5*X(t) + e5(t)
% Y6(t) = gam6*X(t) + e6(t)
% Y7(t) = gam7*X(t) + e7(t)
% Y8(t) = gam8*X(t) + e8(t)
% Y9(t) = gam9*X(t) + e9(t)
% Y10(t) = gam10*X(t) + e10(t)

% where ei(t) ~ N(0,sig2i) for i = 1,2,3,...,10

clc;
clear; 

%% Step 1: DGP
ym = readmatrix("KOR_YC_2021_08_M.xlsx", "Range","b2:k249");

%% Step 2: Principal Component Analysis
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

%% Plot
% Level, Slope, Curvature
Level = Y(:,end);
Slope = Y(:,1) - Y(:,end);
Curvature = (2*Y(:,5) - Y(:,end) - Y(:,1))/2;

figure
subplot(3, 1, 1)
plot([Level PCm(:,1)]);
title('Level with Estimation Plot')

subplot(3, 1, 2)
plot([Slope PCm(:,2)]);
title('Slope with Estimation Plot')

subplot(3, 1, 3)
plot([Curvature PCm(:,3)]);
title('Curvature with Estimation Plot')

