%% Final Exam
% Macro AFNS Model 

clear;
clc;

%% Load Data
YC= readmatrix('KOR_YC_2021_08_Y.xlsx','Range','B2:K249');
Macro = readmatrix('KOR_YC_2021_08_M','Range','B2:C249');

% Raw Data Plot
figure
plot(YC)
title('Raw YC')

figure
plot(Macro)
title('Raw Macro')

% Demean the Raw data
D_YC = demeanc(YC);
D_Macro = detrend(Macro);

% Processing Data Plot
figure
plot(D_YC)
title('Demean YC')

figure
plot(D_Macro)
title('Demean Macro')

% Parameter
tau = [3 6 9 12 18 24 30 36 60 120]';
[T,N] = size(D_YC);
k = 3;
m = 2;
km = k+m;

% Final Data
DF = [D_YC D_Macro];
ym = DF;

%% initial blocking scheme, global
k2 = k^2;
km2 = km^2;
m2 = m^2;

nb1 = k;            % lambda
nb2 = k;            % V
nb3 = k*(k-1)/2;    % Gam
nb4 = 1;            % GQ
nb5 = km2;           % G
nb6 = rows(tau);    % Sig
nb7 = m2; % Omega_mm
nb = [nb1; nb2; nb3;nb4; nb5; nb6; nb7];

nbmh = rows(nb);
nmh = sumc(nb);  % # of parameters in rndblocks

upp = cumsum(nb);
low = [0;upp(1:length(nb)-1)] + 1;

indv = 1:nmh;
indv = indv';

indlambda = indv(low(1):upp(1));
indV = indv(low(2):upp(2));
indGam = indv(low(3):upp(3));
indGQ = indv(low(4):upp(4));
indG  = indv(low(5):upp(5));
indSig = indv(low(6):upp(6));
indOmega_mm = indv(low(7):upp(7));

Spec.k = k;
Spec.tau = tau;
Spec.ym = ym;

Spec.indG = indG;
Spec.indGQ = indGQ;
Spec.indlambda = indlambda;
Spec.indV = indV;
Spec.indGam = indGam;
Spec.indSig = indSig;
Spec.indOmega_mm = indOmega_mm;

psi0 =  zeros(nmh,1);
%% Initial Parameter

% kappa
kappa = -log(0.9); 

% G
G_ = zeros(km2, 1);
G_(1:(km+1):km2) = [0.9, 0.8,kappa*exp(-kappa), 0.7, 0.7];

beta = [1; 1; 0];
GQ = [1 0 0 ; 0 exp(-kappa) kappa*exp(-kappa) ; 0 0 exp(-kappa)];

% lambda
SR = meanc(YC(:,1));   
LR = meanc(YC(:,10));  
m_spread = LR - SR;       
B = makeB(beta, GQ, tau, k);
BL = B(1,:);

lambda = zeros(k,1);
lambda(1,1) = -0.6*m_spread/BL(1,1);
lambda(2,1) = 12*lambda(1);
lambda(3,1) = 20*lambda(1);

% Using PCA

% Data
Y = YC;
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

Subset = PCm( : , 1:3);

%% OLS VAR

[~, Omega_hat1, ~, ~, ~] = OLS_VAR(Subset, 1);

[G_decom,L_decom] = ldl(Omega_hat1); % LGL Decomposition

L_decom(tril(true(size(L_decom)), -1) & L_decom ~= 0)

[~, Omega_hat2, ~, ~, ~] = OLS_VAR(D_Macro, 1);

Omega_mm_hat = Omega_hat2;

%% Starting values at the prior mean
% Beta and Delta
beta = [1; 1; 0];
delta = meanc(ym(:,1));
Spec.beta = beta;
Spec.delta = delta;

psi0(indlambda) = lambda;             % lambda
psi0(indV) = diag(L_decom);
psi0(indGam) =  G_decom(tril(true(size(G_decom)), -1) & G_decom ~= 0);
psi0(indGQ) = kappa;  % GQ
psi0(indG) = G_;
psi0(indSig) = 1*ones(rows(indSig),1);  % Sigma
psi0(indSig) = psi0(indSig);
psi0(indOmega_mm) = vec(Omega_mm_hat);

%% Iteration

load psimx.txt
psi0 = psimx;

%load psi_Con.txt
%psi0 = psi_Con;

% Index
indbj = 1:rows(psi0);
indbj = indbj';

% printi = 1 => See the opimization produdure
% printi = 0 => NOT see the opimization produdure
printi = 1;

% Optimization
[psimx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,psi0,Spec,printi,indbj);
save psimx.txt -ascii psimx;

% Estimates by Deltamethod
thetamx = makeTheta(psimx,Spec);                % Transform psi -> theta
grad = Gradpnew1(@makeTheta,psimx,indbj,Spec);  % Gradient
cov_fnl = grad*Vj*grad';                        % Covariance Matrix
diag_cov = diag(cov_fnl);                       % Variance (diagonal)
stde = sqrt(diag_cov);                          % Standard deviation
t_val = thetamx./stde;                          % t-value
p_val = 2*(1 - cdf('t',abs(t_val),T-k));        % p-value

%결과보기
output_para = [indbj thetamx(indbj) t_val(indbj) stde(indbj) p_val(indbj)];
disp('=======================================');
disp(['Index ', 'Estimates ', ' t value', ' s.e. ', ' p value ']);
disp('------------------------------------------------------------------');
disp(output_para);
disp('------------------------------------------------------------------');
% fm = filtered factors, Pm = condtional variance of factors
% fittedm = fitted values, Residm = residuals

theta = makeTheta(psimx, Spec);
lambda = makeLambda(theta, Spec);
delta = makeDelta(theta, Spec);
G = makeG(theta, Spec);
V = makeV(theta, Spec);
Gam = makeGam(theta, Spec);
Omega = makeOmega(V, Gam, theta,Spec);
L = makeL(Omega);
GQ = makeGQ(theta, Spec);
beta  = makeBeta(theta, Spec);
[a, b, A, B] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, Spec.tau, Spec.k);
Sigma = makeSigma(theta, Spec);
[Fm, Pm, Fittedm, Residm] = KM_filter(a, b, G, Omega, Sigma, ym);

%% Calculate the Term Premium
% Target Term Premium : tau = 60

total_tau = 120;

Gff = G(1:3,1:3);
Gfm = G(1:3,4:5);

Phi_f = inv(L)*(Gff-GQ);
Phi_m = inv(L)*Gfm;
 
Ft = Fm(:,1:3)';
Mt = Fm(:,4:5)';

exr_t = zeros(T,120);

for t = 2:total_tau
    for i = 1:T
        exr_t(i,t) = (-0.5*B(:,t-1)'*L*L'*B(:,t-1))-(B(:,t-1)'*L*lambda)-(B(:,t-1)'*L*Phi_f*Ft(:,i))-(B(:,t-1)'*L*Phi_m*Mt(:,i)); 
    end
end    

Target_Tau = 60;
Term_Premium = 1/60*(sum(exr_t(:, 1:Target_Tau), 2));

%% Result Plot
intvl = 1/12;
startday = 2001+intvl;
endday = 2001+intvl*T;
datat = startday:intvl:endday;
datat = datat';
datat = datat(1:T);

figure
plot(datat, D_YC(:,end), 'b-', datat, Fm(:,1), 'k:', 'Linewidth',2);
legend('Long rates','Level')

figure
plot(datat, D_YC(:,1) - D_YC(:,end), 'b-', datat, Fm(:,2), 'k:' , 'Linewidth',2);
legend('Spread','Slope')

figure
plot(datat, 2*D_YC(:, 5) - (D_YC(:, 1) + D_YC(:, end)), 'b-', datat, Fm(:,3), 'k:', 'Linewidth',2);
legend('Difference between Spreads','Curvature')

figure
plot(datat,D_YC(:,1), 'k:',datat,Fittedm(:,1),'b-', datat,Residm(:,1),'r--', 'Linewidth',2);
legend('Actual','Fitted','Residual')
title('Short rate')

figure
plot(datat,D_YC(:,end),'k:',datat,Fittedm(:,10),'b-', datat,Residm(:,10),'r--', 'Linewidth',2);
legend('Actual','Fitted','Residual')
title('Long rate')

figure
plot(datat, Term_Premium, 'b', 'Linewidth', 2);

