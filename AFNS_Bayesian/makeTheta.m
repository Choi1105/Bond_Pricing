function theta = makeTheta(psi, Spec)

theta = psi;
sf_Sig = Spec.sf_Sig;
sf_Omega = Spec.sf_Omega;
sf_lambda = Spec.sf_lambda;
theta(Spec.ind_Sig) = psi(Spec.ind_Sig)/(sf_Sig^2);
theta(Spec.ind_V) = psi(Spec.ind_V)/sf_Omega;
theta(Spec.ind_lambda) = psi(Spec.ind_lambda)/sf_lambda;
theta(Spec.ind_Omega_MZ) = psi(Spec.ind_Omega_MZ)/(sf_Omega^2);    

end

