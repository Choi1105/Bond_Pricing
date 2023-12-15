function lnL = Kalman_new(ind_B, ind_NB, a, b, G, Omega, Sigma, Ym)

ntau = rows(a);
YCm = Ym(:, 1:ntau);
Macrom = Ym(:, ntau+1:end);

Fm = Gen_Fm_sub(YCm, ind_B, a, b);

T = rows(YCm);
k = length(ind_B);
bbar = b(ind_B, 1:k);
ahat = a(ind_NB);
bhat = b(ind_NB, 1:k);
bbar_inv = inv(bbar);
lnL = 0;
Xm = [Fm, Macrom];
kmz = cols(Xm);

Yhatm = kron(ones(T,1), ahat') + Fm*bhat';
Ehatm = YCm(:, ind_NB) - Yhatm; % Non basis 수익률의 잔차항 
n_NB = length(ind_NB);
Uhatm = Xm(2:end, :) - Xm(1:end-1, :)*G'; % Non basis macro, basis yield의 잔차항 
for t = 2:T
    
    % basis, macro, basis yield의 밀도
    lnpdf_xt = lnpdfmvn(Uhatm(t-1,:)', zeros(kmz,1), Omega);

    % non-basis의 경우 다 독립이기 때문에 개별적으로 계산해서 합
    lnpdf_yhat = sumc(lnpdfn(Ehatm(t, :)', zeros(n_NB, 1), Sigma));
    lnL = lnL + lnpdf_xt + lnpdf_yhat(1,1);
    
end

% 자료의 갯수만큼 jacobian을 더해줘야함
lnL = lnL + (T-1)*log(det(bbar_inv));

end