function [Omega_hat, bbar, fm, phi_hat, B] = make_Omegahat(YCm, Macrom, beta, GQ, ind_basis, tau)

k = rows(beta);

B = makeB(beta, GQ, tau(end), k); % B 만들기, k by max(tau)
B = B';  % max(tau) by k
bbar = B(tau, :)./kron(ones(1, k), tau); % ntau by k

ycm = demeanc(YCm); % 평균이 제거된 수익률 (Y - A(tau))
bbar_basis = bbar(ind_basis, :); % 3 by 3
fm = ycm(:, ind_basis)*inv(bbar_basis)';

% figure
% plot(fm*1200)
xm = demeanc([fm, Macrom]);

y0 = xm(2:end, :); % T-p by k
y_lag = xm(1:end-1, :); % 1기 시차 변수

y_lag2 = y_lag'*y_lag;
phi_hat = y_lag2\(y_lag'*y0); % p*k by k

% Omega_hat
u_hat = y0 - y_lag*phi_hat; % T-p by k
u_hat = demeanc(u_hat);
Omega_hat = u_hat'*u_hat/rows(y0);  % k by k

end