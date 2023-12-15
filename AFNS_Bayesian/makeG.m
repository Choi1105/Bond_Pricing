% Making G matrix
function [G] = makeG(theta,Spec)

ind_G = Spec.ind_G;
kmz = sqrt(rows(ind_G));
G = reshape(theta(Spec.ind_G), kmz, kmz);

end