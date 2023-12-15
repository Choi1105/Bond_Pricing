clear
clc

%% Options
lag = 2;                                      % 자료의 MA 시차 
tau = [3 6 9 12 18 24 30 36 60 120]';         % 만기 
maxac = 50;                                   % 비효율성 계수의 시차
freq = 100;                                   % 중간결과 
n1 = 1000;                                    % 샘플링 크기
n0 = 200;                                     % 번인 크기


%% Data Loading
Macrom_0 = xlsread('KOR_YC_2021_08_M','Macro','B2:C249');
YCm_0 = xlsread('KOR_YC_2021_08_M','YC','B2:K249');
m = cols(Macrom_0);

YCm = YCm_0;
Macrom = Macrom_0;
if lag > 0
    for t = lag+1:length(Macrom_0)
        mt = meanc(Macrom_0(t-lag:t, 1:m));
        Macrom(t, 1:m) = mt';
    end
end

% Macrom = detrend(Macrom);   % 추세제거시 사용
% 1200 나눠주는 이유 = 이론적인 모형과 부합
% Decimal and Monthly 
Macrom = 1*standdc(Macrom(lag+1:end, :))/1200;     % 표준화 or 평균제거
YCm = YCm(lag+1:end, :)/1200;

%% 추정
% 거시변수 설정 
Macrom_M = Macrom(:, 1); % FRP에 영향을 주는 Unspanned Risk 
Macrom_Z = Macrom(:, 2); % FRP에 영향을 주지 않는 Unspanned Risk

% 사전 설정
Spec = make_Spec(YCm, Macrom_M, Macrom_Z, tau, freq, maxac, n0, n1);

k = Spec.k;
theta = Spec.theta;

[lambda, delta, G, Omega, L, GQ, beta] = make_Para(theta, Spec);
Omega_ff = Omega(1:3, 1:3);
[Fm, a, b, A, B, termp] = Gen_Fm(theta, Spec);

Risk_const = zeros(max(tau), 1);
for tauj = 2:max(tau)
    
   B_j =  B(:, 1:tauj-1); % B_full = 3 by 120
   Jensen = 0.5*diag(B_j'*Omega_ff*B_j);
   Risk_const(tauj) = 1200*sumc(-B_j'*lambda - Jensen)/tauj;
   
end

Risk_const = Risk_const(tau);

Omega_ff = Omega_ff*1200*1200;
G(1,1) = 0.9999;
R0 = makeR0(G(1:k, 1:k)', Omega_ff);

ntau = length(tau);
Prior_YCm = zeros(n1, ntau);

for iter = 1:n1
    factor = chol(R0)'*randn(k, 1);
    yc = 1200*a + b*factor;
    Prior_YCm(iter, :) = yc';
end

p = [0.05 0.5 0.95];
Prior_YC = quantile(Prior_YCm,p)';

tau = Spec.tau;
xtick = [1; 12; 36; 60; 84; 120];
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/3])
subplot(1,2,1);
h = plot(tau, meanc(Prior_YCm),'k-', tau, Risk_const, 'b--');
set(gca,'XTick',xtick);
ylim([-0.2 6.5])
xlabel('maturity (month)','FontSize',10)
ylabel('(%)','FontSize',10)
set(h,'LineWidth',2)
legend('yield curve', 'term premium')
title('(a) average yield curve and term premium')

subplot(1,2,2);
h = plot(tau, Prior_YC(:, 1),'b--',Spec.tau, Prior_YC(:, 2),'k-',Spec.tau, Prior_YC(:, 3),'b--');
set(gca,'XTick',xtick);
legend('5%', '50%', '95%')
set(h,'LineWidth',2)
xlabel('maturity (month)','FontSize',10)
ylabel('(%)','FontSize',10)
ylim([-35 35])
title('(b) quantiles of yield curve ')


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3])
h = plot(tau, Risk_const, 'b--');
set(gca,'XTick',xtick);
ylim([-0.2 2])
xlabel('maturity (month)','FontSize',10)
ylabel('(%)','FontSize',10)
set(h,'LineWidth',2)
legend('term premium','location', 'southeast')