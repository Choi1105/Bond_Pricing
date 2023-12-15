

%% Principal Component Analysis
% DGP
% X(t) = mu + phi*X(t-1) + v(t)
% where v(t) ~ N(0,sig2v)

% Y1(t) = gam1*X(t) + e1(t)
% Y2(t) = gam2*X(t) + e2(t)
% Y3(t) = gam3*X(t) + e3(t)
% where ei(t) ~ N(0,sig2i) for i = 1,2,3

clc;
clear; 

%% Step 1: DGP
T = 100; 

% DGP for X
mu = 0.5; 
phi = 0.8;
sig2v = 2;

xm = zeros(T,1);

xm(1) = mu/(1-phi) + randn(1,1)*sqrt(sig2v/(1-phi^2)); 

for t = 2:T
    
    xm(t) = mu + phi*xm(t-1) + randn(1,1)*sqrt(sig2v); 
    
end

% DGP for Yi's

gam1 = -2; 
gam2 = -0.5; 
gam3 = 1;

sig21 = 0.1;
sig22 = 0.2;
sig23 = 0.3;

y1 = gam1*xm + randn(T,1)*sqrt(sig21); 
y2 = gam2*xm + randn(T,1)*sqrt(sig22); 
y3 = gam3*xm + randn(T,1)*sqrt(sig23);

ym = [y1,y2,y3]; 

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

% PCA and Common Factor
figure
plot([PCm(:,1) demeanc(xm)]);