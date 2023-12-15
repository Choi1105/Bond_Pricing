% Unconditional variance of the factors 
function P_ll = makeP_ll(G,omega)

k = rows(G);
k2 = k^2;
G2 = kron(G,G);
eyeG2 = eye(k2)-G2;
omegavec = reshape(omega,k2,1);
P_ll = (eyeG2)\omegavec; 
P_ll = reshape(P_ll,k,k)';
P_ll = (P_ll + P_ll')/2;

end
