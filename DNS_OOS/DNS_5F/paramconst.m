% the parameter constraints
% valid = 0 if invalid
% valid = 1 if valid

function [valid] = paramconst(psi,Sn)

validm = ones(30,1);
if minc(isfinite(psi)) == 0 
    validm(30) = 0;
end

theta = maketheta(psi,Sn);

% Stationarity of the factors
G = reshape(theta(Sn.indG),Sn.k,Sn.k);
eigG = eig(G); 
validm(1) = maxc(abs(eigG)) < 1;

% unconditional variance of latent factors (R0) should be positive definite
Omegam = theta(Sn.indOmega);
Omega = [Omegam(1) Omegam(2) Omegam(3)  Omegam(4)  Omegam(5) ;...
         Omegam(2) Omegam(6) Omegam(7)  Omegam(8)  Omegam(9) ;...
         Omegam(3) Omegam(7) Omegam(10) Omegam(11) Omegam(12);...
         Omegam(4) Omegam(8) Omegam(11) Omegam(13) Omegam(14);...
         Omegam(5) Omegam(9) Omegam(12) Omegam(14) Omegam(15)];
eigOmega = eig(Omega); 
validm(2) = minc(eigOmega) > 0; 

valid = minc(validm); % if any element is equal to zero, invalid

end