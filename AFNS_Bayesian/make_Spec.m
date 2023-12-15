function Spec = make_Spec(YCm, Macrom_M, Macrom_Z, tau, freq, maxac, n0, n1)

Spec.freq = freq;
Spec.maxac = maxac;
Spec.n0 = n0;
Spec.n1 = n1;
Spec.YCm = YCm;

% 수익률 곡선의 평균, 중간결과로 a와 비교하기 위해
Avg_YC = meanc(YCm);
Spec.Avg_YC = Avg_YC;

k = 3;
ntau = length(tau);

m = cols(Macrom_M);
if m == 0
    Macrom_M = [ ];
else
    Macrom_M = demeanc(Macrom_M);
end

z = cols(Macrom_Z);
if z == 0
    Macrom_Z = [ ];
else
    Macrom_Z = demeanc(Macrom_Z);
end

kmz = k + m + z;
Spec.kmz = kmz;
mz = m + z;
Spec.mz = mz;

Spec.Macrom_M = Macrom_M;
Spec.Macrom_Z = Macrom_Z;
Macrom = [Macrom_M, Macrom_Z];
Spec.Macrom = Macrom;

Spec.Ym = [YCm, Macrom];

% SF = scale factor
% 해당 파라미터는 숫자가 너무 작아서 1200을 곱해서 sampling
Spec.m = m;
Spec.z = z;
Spec.tau = tau;
sf_Omega = 1200; 
Spec.sf_Omega = sf_Omega;
sf_Sig = 1200;
Spec.sf_Sig = sf_Sig;
sf_lambda = 1200;
Spec.sf_lambda = sf_lambda;
Spec.tau = tau;
Spec.delta = meanc(YCm(:, 1)); 
% 단기금리의 unconditional mean
% 단기금리는 persistent하기 때문에 I(1)에 가까워
% 추정하게 되면 효율성이 매우 떨어져서 일반적으로 단기금리의 평균으로 고정
% 만약 이 파라미터를 추정하게 되면 다른 파라미터의 효율성이 문제가 발생 
% 수익률 곡선이 I(1)이라고 가정하는 모형도 존재


% Basis Yield 선정에 있어 단기와 중기 설정에 대한 의견은 다름
% 수익률 데이터의 2번째만기를 단기, 6번째만기를 중기, 장기는 10년으로 설정
shortterm = 2;
midterm = 6;

% Basis index
ind_B = [shortterm; midterm; ntau];

% Non-Basis index
if shortterm == 2
    ind_NB0 = 1;
    ind_NB1 = shortterm+1:(midterm-1);
    ind_NB2 = (midterm+1):(ntau-1);
    ind_NB = [ind_NB0; ind_NB1'; ind_NB2'];
else
    ind_NB1 = 2:(midterm-1);
    ind_NB2 = (midterm+1):(ntau-1);
    ind_NB = [ind_NB1'; ind_NB2'];
end
Spec.ind_B = ind_B;
Spec.ind_NB = ind_NB;

% initial blocking scheme, global
kmz2 = kmz^2;
nb1 = k;          % lambda
nb2 = 1;          % GQ
nb3 = k;          % Omega_V
nb4 = k*(k-1)/2;  % Omega_Gam
nb5 = mz^2;       % Omega_MZ
nb6 = kmz2;       % G
nb7 = 1;          %  Sig
% 만기별로 measurement error는 동일하다고 가정
% 차지하는 비중이 작아서 추정결과에 영향이 없음

nb = [nb1; nb2; nb3; nb4; nb5; nb6; nb7];
nmh = sumc(nb);  % # of parameters in rndblocks

upp = cumsum(nb);
low = [0;upp(1:length(nb)-1)] + 1;

indv = 1:nmh;
indv = indv';

ind_lambda = indv(low(1):upp(1));
ind_GQ = indv(low(2):upp(2));
ind_V = indv(low(3):upp(3));
ind_Gam = indv(low(4):upp(4));
ind_Omega_F = [ind_V; ind_Gam];
ind_Omega_MZ = indv(low(5):upp(5));
ind_G  = indv(low(6):upp(6));
ind_Sig  = indv(low(7):upp(7));

Spec.k = k;
Spec.tau = tau;
Spec.kmz = kmz;

Spec.ind_G = ind_G;
Spec.ind_GQ = ind_GQ;
Spec.ind_lambda = ind_lambda;
Spec.ind_V = ind_V;
Spec.ind_Gam = ind_Gam;
Spec.ind_Omega_F = ind_Omega_F;
Spec.ind_Omega_MZ = ind_Omega_MZ;
Spec.ind_Sig = ind_Sig;
Spec.ind_Omega = [ind_Omega_F; ind_Omega_MZ];

beta = [1; 1; 0];
Spec.beta = beta;

% prior for kappa
Spec.kappa_ = -log(0.935); % mean
Spec.kappa_V = 0.00001^2;  % variance

% prior for G
GQ = makeKappa(Spec.kappa_);
G_ = diag([1; 0.935; 0.935; 0.8*ones(m+z,1)]);
G_(1:k, 1:k) = GQ';
G_ = vec(G_);
Gv_ = 0.04*ones(kmz, kmz);
Gv_(1:k, 1:k) = 0.00004*ones(k, k);

if z > 0
    Gv_(1:k,k+m+1:end) = 0.000000000001;
end
Gv_ = vec(Gv_');

Spec.G_ = G_;
Spec.Gv_ = Gv_;

% Make Omega and B by OLS
[Omega_hat, ~, ~, ~, B] = make_Omegahat(YCm, Macrom, beta, GQ, Spec.ind_B, tau);

% prior for lambda
% 평균 수익률 곡선의 장단기차의 약 60%가 기간 프리미엄
% 기간프리미엄은 관측이 안되는데 lambda에 대한 정보가 너무 작음
% 따라서 lambda는 수익률 곡선 모형중에 가장 추정하기 어려움
% 결국 lambda에 대해선 강한 prior를 줘야함
avg_spead = Avg_YC(end) - Avg_YC(1); 
BL = B(end-1, :);
lambda_ = zeros(nb1, 1);
lambda_(1) = -0.6*avg_spead/BL(1,1);
if nb1 > 1
   lambda_(2) = 12*lambda_(1);
elseif nb1 > 2
   lambda_(3) = 20*lambda_(1);
end

Spec.lambda_ = sf_lambda*lambda_; % mean
Spec.lambda_V = (0.001^2)*ones(rows(ind_lambda), 1); % variance

% prior for Omega
% OLS로 추정된 factor에 해당되는 분산-공분산 행렬 
Var_F_hat = diag(Omega_hat(1:k, 1:k));

% m(t), z(t)
Spec.nu0_MZ = 50;
Spec.R0_MZ = inv(diag(diag(Omega_hat(k+1:end, k+1:end)))*(sf_Omega^2))/Spec.nu0_MZ;

% f(t)
Spec.nu0_F = 50;
Spec.R0_F = inv(diag(Var_F_hat)*(sf_Omega^2))/Spec.nu0_F;

% prior for Sigma
Spec.v0 = 50;
Spec.d0 = (0.05^2)*Spec.v0;

% 사전분포가 정규분포인 파라미터들의 사전 평균과 분산
Spec.mu_ = [Spec.lambda_; Spec.kappa_];
Spec.Var_ = [Spec.lambda_V; Spec.kappa_V];

%% 초기값
theta = zeros(nmh, 1);
theta(ind_G) = G_;
theta(ind_lambda) = Spec.lambda_/sf_lambda;
theta(ind_V) = sqrt(Var_F_hat);
theta(ind_GQ) = Spec.kappa_;
theta(ind_Omega_MZ) = vec(diag(diag(Omega_hat(k+1:end, k+1:end))));
theta(ind_Sig) = (0.05^2)/(sf_Sig^2);  % Sigma
Spec.theta = theta;

end

