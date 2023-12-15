function [Sigma] = makeSigma(theta, Spec)

Sigma = theta(Spec.indSig);
k = Spec.k;

Sigma = diag(Sigma);
%Sigma = 0.5*(Sigma + Sigma');

end