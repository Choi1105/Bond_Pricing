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
Rho = theta(Sn.indGamma);
validm(2) = abs(Rho(1)) < 1 ;
validm(3) = abs(Rho(2)) < 1 ;
validm(4) = abs(Rho(3)) < 1 ;

valid = minc(validm); % if any element is equal to zero, invalid

end