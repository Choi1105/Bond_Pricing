% Making L matrix
function L = makeL(Omega)

L = chol(Omega(1:3, 1:3))';

end
 

