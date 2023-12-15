% Likelihood function
function [lnL] = lnlik(psi, Spec)

theta = makeTheta(psi, Spec);
[lambda, delta, G, Omega, ~, GQ, beta, Sigma] = make_Para(theta, Spec);

k = Spec.k;
tau = Spec.tau;
L = cholmod(Omega(1:k,1:k));
[a, b, ~ , ~ , ~] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, tau, k);

Ym = Spec.Ym;
ind_B = Spec.ind_B;
ind_NB = Spec.ind_NB;

% ������ Percent�� ����ϱ�
% scale�� �ʹ� ������ MLE�ϴ� ��쿡 Numerical ���� �߻��ϴ� ��� ���� 
sf_Omega = Spec.sf_Omega;
sf_Sig = Spec.sf_Sig;
Omega = Omega*(sf_Omega^2);
Sigma = Sigma*(sf_Sig^2);
Ym = Ym*1200;
a = a*1200;

lnL = Kalman_new_mex(ind_B, ind_NB, a, b, G, Omega, Sigma, Ym);

end
