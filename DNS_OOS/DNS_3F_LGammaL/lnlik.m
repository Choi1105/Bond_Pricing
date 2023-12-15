%% Likelihood

function [lnL] = lnlik(psi,Sn) 

theta = maketheta(psi,Sn) ; 
[Lam, Sigma, mu, G, L, Gamma] = makePara(theta, Sn);

C = zeros(rows(Lam),1); 
H = Lam;
R = Sigma;

Mu = mu - G*mu; % mu or (demean) mu - G*mu
F = G;
Q = L*Gamma*L; 

ym = Sn.Ym;

lnL = Kalman(C,H,R,Mu,F,Q,ym);

end