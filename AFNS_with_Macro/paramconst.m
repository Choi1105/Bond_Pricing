% the parameter constraints
% valid = 0 if invalid
% valid = 1 if valid
function [valid, validm] = paramconst(psi, Spec)

indV = Spec.indV;
indGam = Spec.indGam;
indGQ = Spec.indGQ;
indSig= Spec.indSig;
indOmega_mm = Spec.indOmega_mm;

validm = ones(30,1);

isfin = isfinite(psi);
validm(1) = minc(isfin) > 0.5;

if validm(1) > 0

    theta = makeTheta(psi, Spec);
    
    validm(2) = minc(theta(indV)) > 0;

    Gam = makeGam(theta, Spec);
    validm(3) = isPositiveDefinite(Gam);
    validm(10) = maxc(abs(theta(indGam))) < 0.95;
    
    validm(4) = maxc(theta(indGQ)) < 1;
    validm(5) = minc(theta(indGQ)) > 0;
   
    G = makeG(theta, Spec);
    eigG = eig(G);
    validm(6) = maxc(abs(eigG)) < 1;
   
    validm(7) = minc(psi(indSig)) > 0;

    O_mm = theta(indOmega_mm);
    validm(8) = O_mm(1) > 0;
    validm(9) = O_mm(4) > 0;
   

    lambda = makeLambda(theta, Spec);
    delta = makeDelta(theta, Spec);
    V = makeV(theta, Spec);
    Omega = makeOmega(V, Gam, theta, Spec);
    L = makeL(Omega);
    GQ = makeGQ(theta, Spec);
    beta  = makeBeta(theta, Spec);

   [a, ~, ~, ~] = makeABbar_HW(delta, G, L, Omega, beta, GQ, lambda, Spec.tau, Spec.k);
   MeanY = meanc(Spec.ym);
   
   dif_yc = a - MeanY(1:10)
   
   if maxc(abs(dif_yc)) < 0.5
       %save psi_Con.txt -ascii psi;
   end
   
   validm(15) = maxc(abs(dif_yc)) < 0.5;

end

valid = minc(validm); % if any element is equal to zero, invalid

end

