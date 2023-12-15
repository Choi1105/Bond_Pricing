function theta = Gen_Sigma1_sub(theta, YCm_NB, Fm, ind_NB, sf, a, b, d0, v0, ind_Sig)


T = rows(Fm);
Fitted = kron(ones(T,1), a(ind_NB)') + Fm*b(ind_NB, :)';
ehatm = YCm_NB - Fitted; % ÀÜÂ÷Ç×

sf2 = sf^2;
d1 = d0;
v1 = v0;
for indtau = 1:rows(ind_NB)
    ehat = sf*ehatm(:,indtau);
    ehat2 = ehat'*ehat;
    d1 = d1 + ehat2;
    v1 = v1 + length(ehat);
end

sig2 = randig(v1/2, d1/2, 1, 1);
theta(ind_Sig) = max(sig2, 0.00000000005^2)/sf2;

end