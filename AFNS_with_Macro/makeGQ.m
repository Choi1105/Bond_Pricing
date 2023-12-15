% Making GQ matrix
function GQ = makeGQ(theta, Spec)

indGQ = Spec.indGQ;
kappa = theta(indGQ);

GQ = makeKappa(kappa);

 
 end