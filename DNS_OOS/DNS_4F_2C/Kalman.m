function [lnL] = Kalman(C,H,R,Mu,F,Q,ym)

T = rows(ym); 
N = cols(ym);

beta_ll = makeBeta_ll(F,Mu); 
P_ll = makeP_ll(F,Q);

lnLm = zeros(T,1); 

for t = 1:T

    beta_tl = Mu + F*beta_ll; % k by 1 
    P_tl=F*P_ll*F'+Q; %k by k
    eta_tl = ym(t,:)' - C - H*beta_tl; % N by 1 
    f_tl = H*P_tl*H' + R; %N by N
    f_tl = (f_tl + f_tl')/2;

    lnLm(t) = lnpdfmvn(eta_tl,zeros(N,1),f_tl);

    Kt = P_tl*H'*invpd(f_tl); % Kalman gain 
    beta_tt = beta_tl + Kt*eta_tl;
    P_tt = P_tl - Kt*H*P_tl;
    beta_ll = beta_tt; P_ll = P_tt;

end

lnL = sumc(lnLm);

if isnan(lnL) == 1; 
    disp('lnL is a NaN'); 
    lnL = -exp(200);
return 

end

end




