% Kalman filter to construct likelihood given sm
function [Fm, a, b, A, B, termp] = Gen_Fm(theta, Spec)

YCm = Spec.YCm;
tau = Spec.tau;
k = 3;
[lambda, delta, G, Omega, L, GQ, beta] = make_Para(theta, Spec);

[a, b, A, B, termp] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, tau, k);

%% Kalman filtering step 
ind_B = Spec.ind_B;
Fm = Gen_Fm_sub_mex(YCm, ind_B, a, b);

end
