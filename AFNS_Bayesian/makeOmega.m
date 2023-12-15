function [Omega, Gam] = makeOmega(theta, Spec)

k = Spec.k;
kmz = Spec.kmz;
mz = kmz - k;

Omega = zeros(kmz, kmz);

V_F = diag(theta(Spec.ind_V));
Gam = triu(0.5*eye(k) + 0.5);
Gam(Gam==0.5) = theta(Spec.ind_Gam); 
Gam = Gam + tril(Gam.',-1);

Omega_F = V_F*Gam*V_F;

Omega(1:k, 1:k) = Omega_F;

Omega_MZ = theta(Spec.ind_Omega_MZ);
Omega(k+1:end, k+1:end) = reshape(Omega_MZ, mz, mz);
Omega = 0.5*(Omega + Omega');

end