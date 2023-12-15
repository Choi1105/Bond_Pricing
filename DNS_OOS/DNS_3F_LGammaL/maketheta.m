function [theta] = maketheta(psi,Sn) 

theta = psi;
theta(Sn.indSig) = exp(psi(Sn.indSig)); 

end