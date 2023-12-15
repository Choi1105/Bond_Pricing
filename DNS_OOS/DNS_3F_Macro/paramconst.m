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
Gm = theta(Sn.indG);
G =  [Gm(1)   Gm(2)   Gm(3)    Gm(4)    Gm(5) ;...
              0         Gm(6)   Gm(7)    Gm(8)    Gm(9) ;...
              0              0        Gm(10) Gm(11) Gm(12)  ;...
              0              0            0          Gm(13)     0  ;...
              0              0            0             0          Gm(14)  ];
eigG = eig(G); 
validm(1) = maxc(abs(eigG)) < 1;


% unconditional variance of latent factors (R0) should be positive definite
Ome = theta(Sn.indOmega);
Omega = [Ome(1) Ome(2) Ome(3) Ome(4) Ome(5) ;...
                    Ome(2) Ome(6) Ome(7) Ome(8) Ome(9) ;...
                    Ome(3) Ome(7) Ome(10) Ome(11) Ome(12) ;...
                    Ome(4) Ome(8) Ome(11) Ome(13) Ome(14) ;...
                    Ome(5) Ome(9) Ome(12) Ome(14) Ome(15) ];
eigOmega = eig(Omega); 
validm(2) = minc(eigOmega) > 0; 

valid = minc(validm); % if any element is equal to zero, invalid

end