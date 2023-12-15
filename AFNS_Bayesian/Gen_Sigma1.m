function [theta, ehatm] = Gen_Sigma1(theta, Spec)

[Fm, a, b] = Gen_Fm(theta, Spec);

YCm = Spec.YCm;
ind_Sig = Spec.ind_Sig;
d0 = Spec.d0;
v0 = Spec.v0;
sf = Spec.sf_Sig;
ind_NB = Spec.ind_NB;
YCm_NB = YCm(:, ind_NB);

theta = Gen_Sigma1_sub_mex(theta, YCm_NB, Fm, ind_NB, sf, a, b, d0, v0, ind_Sig);

end