% Making GQ matrix
function GQ = makeGQ(theta, Spec)

ind_GQ = Spec.ind_GQ;
kappa = theta(ind_GQ);

GQ = makeKappa(kappa);

 
 end