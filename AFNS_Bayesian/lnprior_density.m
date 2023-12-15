% compute the log prior density
function lnpriord = lnprior_density(theta, Spec)

% lambda
sf_lambda = Spec.sf_lambda;
 lnpriord = sumc(lnpdfn(sf_lambda*theta(Spec.ind_lambda), Spec.lambda_, Spec.lambda_V));
 
 % Omega_LL
 k = Spec.k;
 sf_Omega = Spec.sf_Omega;
 sf_Sig = Spec.sf_Sig;
 Omega = makeOmega(theta, Spec);
 Omega = (sf_Omega^2)*Omega;
 Omega_inv_F = invpd(Omega(1:k,1:k));
 lnprior_Omega_F = lnpdfWishart(Omega_inv_F, Spec.R0_F, Spec.nu0_F);
 lnpriord = lnpriord + lnprior_Omega_F;
 
 %% kappa
 lnpriord = lnpriord + lnpdfn(theta(Spec.ind_GQ), Spec.kappa_, Spec.kappa_V);

 %% G
 vecG = theta(Spec.ind_G);
 G_ = Spec.G_;
 Gv_ = Spec.Gv_;
 vecG = vecG(abs(Gv_) > 0.00001);
 G_ = G_(abs(Gv_) > 0.00001);
 Gv_ = Gv_(abs(Gv_) > 0.00001);
 lnpriord_G = sumc(lnpdfn(vecG, G_, Gv_));
 lnpriord = lnpriord + lnpriord_G;
 
 %% Omega_MZ
 Omega_inv_MZ = invpd(Omega(k+1:end,k+1:end));
 lnprior_Omega_MZ = lnpdfWishart(Omega_inv_MZ, Spec.R0_MZ, Spec.nu0_MZ);
 lnpriord = lnpriord + lnprior_Omega_MZ;
 
  %% Sigma
 lnprior_Sig = lnpdfig((sf_Sig^2)*theta(Spec.ind_Sig), Spec.v0/2, Spec.d0/2);
 lnpriord = lnpriord + lnprior_Sig;
 
 
 if isnan(lnpriord) == 1
    lnpriord = -exp(20);
 end

end