%% Likelihood

function [lnL] = lnlik(psi,Sn) 

theta = maketheta(psi,Sn) ; 
[Lam, Sigma, mu, G, Omega] = makePara(theta, Sn);

C = zeros(rows(Lam),1); 
H = Lam;
R = Sigma;

Mu = mu - G*mu; % mu or (demean) mu - G*mu
F = G;
Q = Omega; 

ym = Sn.Ym;

lnL = Kalman_mex(C,H,R,Mu,F,Q,ym);

end