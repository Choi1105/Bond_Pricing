% compute the log prior density
function lnpriord = lnprior(psi, Spec)

 lnpriord = sumc(lnpdfn(psi(Spec.ind_lambda), Spec.lambda_, Spec.lambda_V));
 
 k = Spec.k;
 Omega = makeOmega(psi, Spec);
 Omega_inv = invpd(Omega(1:k,1:k));
 lnprior_Omega_F = lnpdfWishart(Omega_inv, Spec.R0_F, Spec.nu0_F);
 lnpriord = lnpriord + lnprior_Omega_F;
 
 lnpriord = lnpriord + lnpdfn(psi(Spec.ind_GQ), Spec.kappa_, Spec.kappa_V);

 if isnan(lnpriord) == 1
    lnpriord = -exp(20);
 end

end