function theta = makeTheta(psi, Spec)

theta = psi;
theta(Spec.indV) = psi(Spec.indV);
theta(Spec.indSig) = psi(Spec.indSig);

end