% Kalman filter

function [Beta_ttm, P_ttm, Fittedm, Residm] = KM_filter(C,H,R,Mu,F,Q,ym)

T = rows(ym); 
N = cols(ym);
k = rows(F);

beta_ll = makeBeta_ll(F,Mu); 
P_ll = makeP_ll(F,Q);

Beta_ttm = zeros(T,k); 
P_ttm = zeros(T,k); 
Fittedm = zeros(T,N); 
Residm = zeros(T,N);

for t = 1:T
    
    beta_tl = Mu + F*beta_ll; % k by 1 % 
    P_tl=F*P_ll*F'+Q; %k by k%
    y_tl = C + H*beta_tl;
    eta_tl = ym(t,:)' - y_tl; % N by 1 %
    f_tl=H*P_tl*H'+R; %N by N%
    f_tl = (f_tl + f_tl')/2;

    Kt = P_tl*H'*invpd(f_tl); 
    beta_tt = beta_tl + Kt*eta_tl; 
    P_tt = P_tl - Kt*H*P_tl;

    Beta_ttm(t,:) = beta_tt'; 
    P_ttm(t,:) = diag(P_tt)'; 
    Pred_Ym = C + H*beta_tt;
    Fittedm(t,:) = Pred_Ym';
    Residm(t,:) = ym(t,:)' - Pred_Ym;

    beta_ll = beta_tt; 
    P_ll = P_tt;

end

end