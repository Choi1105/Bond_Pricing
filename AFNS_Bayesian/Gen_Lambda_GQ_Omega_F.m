function [theta, lnlik0, lnprior0, lnpost0, accept] = Gen_Lambda_GQ_Omega_F(theta,lnpost, lnlik,lnprior,paramconst,Spec)

index = Spec.ind_MH;
psi = makePsi(theta, Spec);
nu = 20; % degree of freedom for multivariate t proposal density

% MH step -- note that only theta is updated in the MH step
[psi, accept, lnlik0, lnprior0, lnpost0] = mhstep(lnpost, lnlik,lnprior,paramconst,psi,index,nu, Spec);
theta = makeTheta(psi, Spec);

end



