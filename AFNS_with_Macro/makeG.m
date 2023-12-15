% Making G matrix
function [G] = makeG(theta,Spec)

indG = Spec.indG;
km = sqrt(rows(indG));
G = reshape(theta(Spec.indG), km, km);

end