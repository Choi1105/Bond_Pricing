% Unconditional mean of the factors
function beta_ll = makeBeta_ll(F,mu)

k = rows(F);
eyeF = eye(k)-F; 
beta_ll = eyeF\mu;

end