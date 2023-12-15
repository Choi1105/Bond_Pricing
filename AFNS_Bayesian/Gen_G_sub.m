function [G, vecG] = Gen_G_sub(Y, Omega_inv, G_, Gv_, G0)

Y = demeanc(Y);
[Y0, YLm] = makeYX(Y);

X = YLm; % 설명변수, 3차원
XX = X(:, :, 1)'*Omega_inv*X(:, :, 1);
XY = X(:, :, 1)'*Omega_inv*Y0(1, :)';
T0 = rows(Y0); % = T-p
k = cols(Y0);
k2 = k^2;

for t = 2:T0
    Xt = X(:, :, t);
    XX = XX + Xt'*Omega_inv*Xt;
    XY = XY + Xt'*Omega_inv*Y0(t,:)';
end

precb_ = diag(ones(k2, 1)./Gv_);
B1_inv = precb_ + XX;
B1_inv = 0.5*(B1_inv + B1_inv');
B1 = inv(B1_inv);
B1 = 0.5*(B1 + B1');
A = XY + precb_*G_; % b_ = B0
BA = B1*A; % full conditional mean

[Chol_B1, is_error] = chol(B1);
if is_error == 0
    vecG = BA + Chol_B1'*randn(k2,1); % beta sampling 하기
    % F 행렬만들기
    G = reshape(vecG, k, k);  % p*k by k
    G = G';
    vecG = vec(G);
    % 안정성 확인하기 크게 중요하지 않음
    eigG = eig(G); % eigenvlaue 계산
    if maxc(abs(eigG)) >= 1
        G = G0;
        vecG = vec(G);
    end
else
    G = G0;
    vecG = vec(G0);
end

end