function [theta, G, vecG] = Gen_G(Y, theta, G_, Gv_, Spec)

G0 = makeG(theta, Spec);
Omega = makeOmega(theta, Spec);
Omega_inv = invpd(Omega);
Omega_inv = 0.5*(Omega_inv + Omega_inv');

[G, vecG] = Gen_G_sub_mex(Y, Omega_inv, G_, Gv_, G0);
k = Spec.k;
m = Spec.m;
z = Spec.z;

if z > 0
  G(1:k, k+m+1:end) = 0;
end

theta(Spec.ind_G) = vec(G);

end