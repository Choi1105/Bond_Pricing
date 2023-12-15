% the parameter constraints
% valid = 0 if invalid
% valid = 1 if valid
function [valid, validm] = paramconst(psi, Spec)

theta = makeTheta(psi, Spec);
ind_GQ = Spec.ind_GQ;
V = psi(Spec.ind_V);
Gam = makeGam(theta, Spec);

validm = ones(30,1);

isfin = isfinite(psi);
validm(1) = minc(isfin) > 0.5;

if validm(1) > 0
    
    validm(7) = minc(psi(ind_GQ)) > 0;
    
    validm(2) = minc(V) > 0.0001;
    validm(3) = isPositiveDefinite(Gam);
    if  validm(3) == 1
        k = Spec.k;
        tau = Spec.tau; 
        [lambda, delta, G, Omega, ~, GQ, beta] = make_Para(theta, Spec);
        L = cholmod(Omega(1:k,1:k));
        [a, ~, ~ , ~ , ~] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, tau, k);
        
        dif_yc = 1200*(a - Spec.Avg_YC);
        %validm(8) = maxc(abs(dif_yc)) < 0.5;
    end
end

valid = minc(validm); % if any element is equal to zero, invalid

end
