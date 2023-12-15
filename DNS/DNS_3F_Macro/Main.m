%% Dynamic Nelson Siegel Model
% 3 latent factor with Macro Variable (L.S.C)

% M.E
% [   y(t)  ] =  [Lam 0    ]  * [ f(t)  ]+  [e(t)]  ,e(t)~iidN(0,Sig)
% [  m(t) ] =  [  0    Im  ]  * [m(t) ] + [  0  ]
% T.E
% f(t) = mu + G*(f(t) - mu) + v(t),v(t)~iidN(0,Omega)

% Lambda Matrix 
% Lam = ones(ntau, 3);
% Lam(:,2) = (ones(ntau,1)-exp(-tau*lambda))./(tau*lambda);
% Lam(:,3) = ((ones(ntau,1)-exp(-tau*lambda))./(tau*lambda)) - exp(-tau*lambda);

% 수정할 .m file  ( makePara.m, paramconst.m, lnlik.m(Kalman_mex -> Kalman으로 수정) )
clear;
clc;
%% Data
% Load Data
YData = readmatrix('Yield.xlsx','Range','B2:K249');
MData = readmatrix('Macro.xlsx','Range','B2:C249');

Data = [YData, MData];

tau = [3 6 9 12 18 24 30 36 60 120]';
[T,N] = size(Data);
Ym = Data;
k = 5; % Factor 개수에 맞게 수정

% Decay Parameter
lambda = 0.0609; % 0.0498, 0.0609, 0.0747

% initial blocking scheme, global
nb1 = N;    % Sigma
nb2 = k;    % Mu
nb3 = (k^2)-11;  % G
nb4 = k*(k+1)/2; % Omega

nb = [nb1;nb2;nb3;nb4];

upp = cumsum(nb);
low = [0;upp(1:rows(nb)-1)] + 1;
indv = 1:sumc(nb);
indv = indv';
indSig = indv(low(1):upp(1));
indMu = indv(low(2):upp(2));
indG = indv(low(3):upp(3));
indOmega = indv(low(4):upp(4));

% initials
sigma0 = 0.01*ones(N,1);
MU0 = 0.01*ones(k,1);
G0 = [0.1;0;0;0;0;0.1;0;0;0;0.1;0;0;0.1;0];
vecG0 = vec(G0);
Omgega0 = [0.1;0;0;0;0;0.1;0;0;0;0.1;0;0;0.1;0;0.1]; % Omega 개수에 맞게 수정
psi0 = [log(sigma0);MU0;vecG0;Omgega0];

%load 'psimx.txt' 
%psi0 = psimx; 

% Structure variables
Sn.indSig = indSig;
Sn.indMu = indMu;
Sn.indG = indG;
Sn.indOmega = indOmega;
Sn.Ym = Ym;
Sn.lambda = lambda;
Sn.k = k;
Sn.tau = tau;

% printi = 1 => See the opimization produdure
% printi = 0 => NOT see the opimization produdure
printi = 1;

% Optimization by Block
% n0 = 10;
% n1 = 50;
% lnLm = zeros(n1,1);
% 
% for iter = 1:n1
%     
%     for blockj = 1:rows(upp)
%         
%         indbj = indv(low(blockj):upp(blockj));
%         [psimx, fmax, ~, ~] = SA_Newton(@lnlik,@paramconst,psi0,Sn,printi,indbj);
%         psi0(indbj) = psimx(indbj);
%         
%     end
%     
%     lnLm(iter) = fmax;
%     
%     if iter > n0
%         if lnLm(iter) - lnLm(iter-1) < exp(-3.5)
%             iter = n1 + 1;
%         end 
%     end
%     
% end

% Optimization
[psimx, fmax,Vj, Vinv] = SA_Newton(@lnlik,@paramconst,psi0,Sn,printi,indv);
save psimx.txt -ascii psimx;

% Estimates by Deltamethod
thetamx = maketheta(psimx,Sn);                  % Transform psi -> theta
grad = Gradpnew1(@maketheta,psimx,indv,Sn);    % Gradient
cov_fnl = grad*Vj*grad';                        % Covariance Matrix
diag_cov = diag(cov_fnl);                       % Variance (diagonal)
stde = sqrt(diag_cov);                          % Standard deviation
t_val = thetamx./stde;                          % t-value
p_val = 2*(1 - cdf('t',abs(t_val),T-k));        % p-value

%결과보기
output_para = [indv thetamx(indv) t_val(indv) stde(indv) p_val(indv)];
disp('=======================================');
disp(['Index ', 'Estimates ', ' t value', ' s.e. ', ' p value ']);
disp('------------------------------------------------------------------');
disp(output_para);
disp('------------------------------------------------------------------');
% fm = filtered factors, Pm = condtional variance of factors
% fittedm = fitted values, Residm = residuals

theta = maketheta(psimx,Sn) ;
[Lam, Sigma, mu, G, Omega] = makePara(theta, Sn);
C = zeros(rows(Lam),1);
H = Lam;
R = Sigma;
Mu = mu - G*mu; % mu or (demean) mu - G*mu
F = G;
Q = Omega;
[Fm, Pm, Fittedm, Residm] = KM_filter(C,H,R,Mu,F,Q,Ym);

intvl = 1/12;
startday = 2001+intvl;
endday = 2001+intvl*T;
datat = startday:intvl:endday;
datat = datat';
datat = datat(1:T);


%% Figure
Level = Ym(:,end);
Slope = Ym(:,1) - Ym(:,end);
Curvature = (2*Ym(:,5) - Ym(:,end) - Ym(:,1))/2;

figure
plot(datat, Level, 'b-', datat, Fm(:,1), 'k:', 'Linewidth',2);
legend('Long rates','Level')

figure
plot(datat, Slope, 'b-', datat, Fm(:,2), 'k:' , 'Linewidth',2);
legend('Spread','Slope')

figure
plot(datat, Curvature, 'b-', datat, Fm(:,3), 'k:' , 'Linewidth',2);
legend('Spread','Curvature')


figure
plot(datat,Ym(:,1), 'k:',datat,Fittedm(:,1),'b-', datat,Residm(:,1),'r--', 'Linewidth',2);
legend('Actual','Fitted','Residual')
title('Short rate')

figure
plot(datat,Ym(:,end),'k:',datat,Fittedm(:,end),'b-', datat,Residm(:,end),'r--', 'Linewidth',2);
legend('Actual','Fitted','Residual')
title('Long rate')

%% RMSE
numRowsToExtract = 12; % 추출할 행의 개수
last12RowsOfEachColumn = Residm(end - numRowsToExtract + 1:end, :);

columnNames = {'3', '6', '9', '12', '18', '24', '30', '36', '60', '120'};
rmseValues = zeros(1, size(last12RowsOfEachColumn, 2));

for col = 1:size(last12RowsOfEachColumn, 2)

    predicted = zeros(12, 1); 
    actual = last12RowsOfEachColumn(:, col); 
    
    % RMSE 계산
    error = predicted - actual;
    squared_error = error.^2;
    mse = mean(squared_error);
    rmseValues(col) = sqrt(mse); 
end

for col = 1:length(columnNames)
    fprintf('RMSE for colname %s:     %.2f\n', columnNames{col}, rmseValues(col));
end