% Parameters in S-S form
function [Lam,Sigma,mu,G,Omega] = makePara(theta, Sn) 

tau = Sn.tau;
lambda = Sn.lambda;

% M.E.
ntau = rows(Sn.tau);
Lam = ones(ntau, 4);
Lam(:,2) = (ones(ntau,1)-exp(-tau*lambda))./(tau*lambda);
Lam(:,3) = (ones(ntau,1)-exp(-tau*3*lambda))./(tau*3*lambda);
Lam(:,4) = ((ones(ntau,1)-exp(-tau*lambda))./(tau*lambda)) - exp(-tau*lambda);

Sigma = diag(theta(Sn.indSig));

% T.E.
mu = theta(Sn.indMu); 
G = reshape(theta(Sn.indG),Sn.k,Sn.k);
Omegam = theta(Sn.indOmega);
Omega = [Omegam(1) Omegam(2) Omegam(3) Omegam(4);...
         Omegam(2) Omegam(5) Omegam(6) Omegam(7);...
         Omegam(3) Omegam(6) Omegam(8) Omegam(9);...
         Omegam(4) Omegam(7) Omegam(9) Omegam(10)];
end