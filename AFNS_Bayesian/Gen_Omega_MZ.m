% Standard

function [theta] = Gen_Omega_MZ(theta, Spec)

Fm = Gen_Fm(theta, Spec);
Xm = [Fm, Spec.Macrom];

T = length(Xm);
sf_Omega = Spec.sf_Omega;
Xt = Xm(2:end, :);
XL = Xm(1:end-1, :);
G = makeG(theta, Spec);
ehat = sf_Omega*(Xt - XL*G'); % T by k

k = Spec.k;
ehat = ehat(:,k+1:end);
ehat2m = ehat'*ehat;
R1_inv = ehat2m + invpd(Spec.R0_MZ);
R1 = invpd(R1_inv);
nu1 = T + Spec.nu0_MZ;

Omega_inv_hat = randwishart(R1, nu1);
Omega_MZ = invpd(Omega_inv_hat);
theta(Spec.ind_Omega_MZ) = vec(Omega_MZ)/(sf_Omega^2);

end