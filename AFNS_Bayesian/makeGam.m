% Making Gamma matrix
function [Gam] = makeGam(theta, Spec)

Gam = triu(0.5*eye(Spec.k) + 0.5);
Gam(Gam==0.5) = theta(Spec.ind_Gam); 
Gam = Gam + tril(Gam.',-1);

end