function [lambda, delta, G, Omega, L, GQ, beta, Sigma] = make_Para(theta, Spec)

lambda = makeLambda(theta, Spec);
delta = makeDelta(theta, Spec);
G = makeG(theta, Spec);
Omega = makeOmega(theta, Spec);
L = makeL(Omega);
GQ = makeGQ(theta, Spec);
beta  = Spec.beta;
Sigma = makeSigma(theta, Spec);

end