function [phi_hat, Omega_hat, F, lnL, BIC] = OLS_VAR(Y, p)

[T, k] = size(Y);
y0 = Y(p+1:T,:); % T-p by k

y_lag = [ ];

for j = 1:p
  y_lag = [y_lag, Y(p+1-j:T-j,:)]; 
end

y_lag2 = y_lag'*y_lag;
phi_hat = y_lag2\(y_lag'*y0); % p*k by k
 
F = [phi_hat'; eye((p-1)*k), zeros(k*(p-1),k)];  % p*k by p*k

T = rows(y0);

% Omega_hat 
u_hat = y0 - y_lag*phi_hat; % T-p by k 
Omega_hat = u_hat'*u_hat/(T-p*k);  % k by k

lnL = 0;
T = rows(u_hat);
    for t = 1:T
        lnL = lnL + lnpdfmvn(u_hat(t, :)', zeros(k, 1), Omega_hat);
    end
    
num_para = p*k^2 + k*(1+k)/2;    
BIC = -2*lnL + num_para*log(T);    
end