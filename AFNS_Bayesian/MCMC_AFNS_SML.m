function [theta_st, post_param, post_summary, lnpostm, lnlikm, Post_FM, Spec] = MCMC_AFNS_SML(Spec)

% Structural Parameters
n0 = Spec.n0;
n1 = Spec.n1;
maxac = Spec.maxac;
freq = Spec.freq;
k = Spec.k;

ind_G = Spec.ind_G;
ind_GQ = Spec.ind_GQ;
ind_lambda = Spec.ind_lambda;
ind_Omega_F = Spec.ind_Omega_F;
ind_Omega_MZ = Spec.ind_Omega_MZ;
ind_Sig = Spec.ind_Sig;

% 자료 변환
Macrom = Spec.Macrom;
T = rows(Macrom);

% MH로 샘플링할 파라미터들
ind_MH = [ind_Omega_F; ind_GQ; ind_lambda];

Spec.ind_MH = ind_MH;
YCm = Spec.YCm;
avg_YC = 1200*meanc(YCm);

%% Starting values at the prior mean
theta = Spec.theta;
theta_st = theta;
lnL = lnlik(theta, Spec);  % log likelihood evaluated at the starting value
lnpriord = lnprior(theta, Spec);
lnpost_st = lnL + lnpriord;

counter = zeros(rows(theta), 1);
post_param = zeros(n0+n1, length(theta));

lnpostm = zeros(n0 + n1, 1);
lnlikm = zeros(n0 + n1, 1);
Post_FM = zeros(n1, T, k);

for iter = 1:(n0+n1)

    %% Omega_MZ 샘플링
    theta = Gen_Omega_MZ(theta, Spec);
    counter(ind_Omega_MZ) = counter(ind_Omega_MZ) +  1;


    %% G Sampling
    [Fm, a] = Gen_Fm(theta, Spec);
    Xm = [Fm, Macrom];
    [theta, G] = Gen_G(Xm, theta, Spec.G_, Spec.Gv_, Spec);
    counter(ind_G) = counter(ind_G) +  1;


    %% Omega_F, kappa, and lambda sampling
    [theta, lnlik0, lnprior0, lnpost0, accept] = Gen_Lambda_GQ_Omega_F(theta,@lnpost, @lnlik,@lnprior,@paramconst,Spec);
    counter(ind_MH) = counter(ind_MH) +  accept;


    %% Sigma 샘플링
    theta = Gen_Sigma1(theta, Spec);
    counter(ind_Sig) = counter(ind_Sig) +  1;

    lnprior0 = lnprior_density(theta, Spec);
    lnpost0 = lnlik0 + lnprior0;
    lnpostm(iter) = lnpost0;
    lnlikm(iter) = lnlik0;
    if lnpost0 > lnpost_st
        theta_st = theta;
        lnpost_st = lnpost0;
    end

    %% 저장 및 디스플레이
    psi = makePsi(theta, Spec);
    post_param(iter, :) = psi';

    if iter > n0
        for j = 1:k
            Post_FM(iter - n0, :, j) = Fm(:, j)';
        end
    end

    if isequal(floor(iter/freq),iter/freq) == 1  % At every 'freq'th iteration,

        clc
        disp('-------------------------------------')
        disp(['current iteration = ', num2str(iter)]);
        disp('-------------------------------------')


        GQ = makeGQ(theta, Spec);
        [~, Gam] = makeOmega(theta, Spec);
        disp('G_ff      GQ')
        disp([G(1:3, 1:3), GQ]);
        disp('G_ff - GQ')
        disp(G(1:3, 1:3)-GQ);
        disp('Gam')
        disp(Gam);

        psi = makePsi(theta, Spec);
        Omega = makeOmega(psi, Spec);
        Omega_MZ = Omega(k+1:end, k+1:end);
        disp(['Vol(MZ) ', num2str(sqrt(diag(Omega_MZ)'))]);
        disp(['diag(Gmm) ', num2str(diag(G(4:end, 4:end))')]);
        disp(['Vol(F) ', num2str(psi(Spec.ind_V)')]);
        disp(['1200*lambda and kappa ', num2str(psi(Spec.ind_lambda)'),'    ', num2str(psi(Spec.ind_GQ))]);
        disp(['sq. Sigma ', num2str(sqrt(psi(Spec.ind_Sig)'))]);

        [Fm, a, ~, ~, ~, termp] = Gen_Fm(theta, Spec);
        
        % 이 둘은 비슷해야함. 달라지면, 모형, 코드, prior문제 일 가능성
        disp(['1200*a = ', num2str(1200*a(Spec.ind_B)')]);
        disp(['ave. YC = ', num2str(avg_YC(Spec.ind_B)')]);

        disp(['ave. TP = ', num2str(termp*1200)]);
        disp(['average factors = ', num2str(1200*meanc(Fm)')]);
        disp(['Acceptance rate = ', num2str(100*counter(ind_GQ)/iter), ' %']);
        disp(['lnlik = ', num2str(lnlik0)]);
        disp(['lnprior = ', num2str(lnprior0)]);
        disp(['lnpost = ', num2str(lnpost0)]);
        disp(['lnpost_st = ', num2str(lnpost_st)]);
        if iter > 2*maxac
            if iter < n0 + 2*maxac
                inef = ineff(post_param(1:iter, :),maxac);
            else
                inef = ineff(post_param(n0+1:iter, :),maxac);
            end
            disp(['ineff. lambda = ', num2str(inef(Spec.ind_lambda)')]);
            disp(['ineff. Omega_F = ', num2str(inef(Spec.ind_Omega_F)')]);
        end
    end
end

accept_rate = counter/(n0+n1);
post_param = post_param(n0+1:end, :);
postmom_para = MHout(post_param, 0.05, 50, 0);

post_summary = [postmom_para, accept_rate*100];

% save theta_st theta_st
lnpostm = lnpostm(n0+1:end);
lnlikm = lnlikm(n0+1:end);

end

