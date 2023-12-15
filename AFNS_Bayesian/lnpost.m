% Likelihood function  
function lnposteriord = lnpost(psi, Spec)

    lnL = lnlik(psi, Spec);
    lnpriord = lnprior(psi, Spec);
    lnposteriord = lnL + lnpriord;
  
end
