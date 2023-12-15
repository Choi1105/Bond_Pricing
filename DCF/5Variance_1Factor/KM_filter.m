%% Kalman Filter
% Wild guess b/c random walk component
% input: C,H,R,Mu,F,Q (ss paramters), ym(data)
% output: lnL (log likelihood), Beta_ttm, P_tt (Filtered values)

function [lnL, Beta_ttm, P_ttm]= KM_filter(C,H,R,Mu,F,Q,ym)

% Data infomation
T = rows(ym); 
N = cols(ym); 
k = rows(F); 

% Initial value for Wild guess
beta_ll = makeBeta_ll(F,Mu); 
P_ll = makeP_ll(F,Q);

% Pre-allocation
lnLm = zeros(T,1); 
Beta_ttm = zeros(T,k);
P_ttm = zeros(T,k); 

% Kalman Filter, refer to Kim and Nelson (1999)
for t = 1:T
    
    % Conditional mean/var for transition eq.
    beta_tl = Mu + F*beta_ll; 
    P_tl = F*P_ll*F'+ Q; 
    
    % Conditional mean/var for measurement eq.
    eta_tl = ym(t,:)' - C - H*beta_tl; 
    f_tl = H*P_tl*H'+ R; 
    f_tl = (f_tl + f_tl')/2; 
    
    % Compute likelihood values
    lnLm(t) = lnpdfmvn(eta_tl,zeros(N,1),f_tl);

    % Kalmain gain = weight on new information
    Kt = P_tl*H'*invpd(f_tl); 
   
    % Updating
    beta_tt = beta_tl + Kt*eta_tl;
    P_tt = P_tl - Kt*H*P_tl; 
    beta_ll = beta_tt;
    P_ll = P_tt;
    
    % Save b(t|t), P(t|t)
    Beta_ttm(t,:) = beta_tt';
    P_ttm(t,:) = diag(P_tt)';
   
end

% Sum of log likelihood
lnL = sumc(lnLm);

% If "NAN", skip!!!
if isnan(lnL) == 1 
    disp('lnL is a NaN'); 
    lnL = -exp(200);  
end

end