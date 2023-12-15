function Omega = makeOmega(V,Gam, theta, Spec)

Omega_mm = theta(Spec.indOmega_mm);
Omega_mmk = reshape(Omega_mm,2,2);

Omega = zeros(5:5);

Omega(1:3,1:3) = V*Gam*V;
Omega(4:5,4:5) = Omega_mmk;
Omega = 0.5*(Omega + Omega');
end