% Likelihood function
function [lnL] = lnlik(psi, Spec)

theta = makeTheta(psi, Spec);
lambda = makeLambda(theta, Spec);
delta = makeDelta(theta, Spec);
G = makeG(theta, Spec);
V = makeV(theta, Spec);
Gam = makeGam(theta, Spec);
Omega = makeOmega(V, Gam, theta,Spec);
L = makeL(Omega);
GQ = makeGQ(theta, Spec);
beta  = makeBeta(theta, Spec);

[a, b, ~, ~] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, Spec.tau, Spec.k);
Sigma = makeSigma(theta, Spec);
ym = Spec.ym;
 
lnL = Kalman(a, b, G, Omega, Sigma, ym);
    % if error, please use Kalman or convert Kalman to Kalman_mex yourself
    % using the C compiler.

end
