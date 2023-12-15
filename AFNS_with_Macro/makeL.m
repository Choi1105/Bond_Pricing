% Making L matrix
function L = makeL(Omega)

L = cholmod(Omega(1:3, 1:3))';

end
 

