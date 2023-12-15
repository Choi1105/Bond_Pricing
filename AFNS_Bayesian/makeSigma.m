function [Sigma] = makeSigma(theta, Spec)

ind_NB = Spec.ind_NB;
Sigma = theta(Spec.ind_Sig)*ones(length(ind_NB), 1);

end