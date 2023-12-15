%% Kalman Smoothing
% input: C,H,R,Mu,F,Q (ss paramters), ym(data)
% output: Beta_ttm, P_tt (Smoothed values)

function [Beta_tTm, P_tTm]= KM_smooth(C,H,R,Mu,F,Q,ym)

% Data infomation
T = rows(ym);
k = rows(Mu);

% Initial value for Wild guess
beta_ll = makeBeta_ll(F,Mu); 
P_ll = makeP_ll(F,Q);

%% Step 1: Forward (Kalman Filter)
% Exactly same as KM_filter.m
Beta_ttm = zeros(T,k);
Beta_tlm = zeros(T,k);
P_ttm = zeros(k,k,T);
P_tlm = zeros(k,k,T);

for t = 1:T
    
    beta_tl = Mu + F*beta_ll;
    Beta_tlm(t,:) = beta_tl';
    
    P_tl = F*P_ll*F'+ Q;
    P_tlm(:,:,t) = P_tl;
    
    eta_tl = ym(t,:)' - C - H*beta_tl;
    f_tl = H*P_tl*H'+ R;
    f_tl = (f_tl + f_tl')/2;
    
    Kt = P_tl*H'*invpd(f_tl);
    beta_tt = beta_tl + Kt*eta_tl;
    P_tt = P_tl - Kt*H*P_tl;
    
    Beta_ttm(t,:) = beta_tt';
    P_ttm(:,:,t) = P_tt;
    
    beta_ll = beta_tt;
    P_ll = P_tt;
    
end

%% Step 2: Backward (Smoothing equation)
% (1) beta_t|T = beta_t|t + P_t|t*F'*P_t+1|t*(beta_t+1|T - beta_t+1|t)
% (2) P_t|T = P_t|t + P_t|t*F'*P_t+1|t*(P_t+1|T - P_t+1|t)*P_t+1|t*F*P_t|t

% Pre-allocation
Beta_tTm = zeros(T,k);
P_tTm = zeros(T,k);

% Compute at T
beta_tt = Beta_ttm(T,:)';
beta_tT = beta_tt;
Beta_tTm(T,:) = beta_tT';
beta_t1T = beta_tt;

P_tt = P_ttm(:,:,T);
P_tT = P_tt;
P_tTm(T,:) = diag(P_tT)';
P_t1T = P_tT;

% Compute at T-1 to 1
t = T - 1;
while t >= 1
    
    % Given Beta_t|t and P_t|t
    beta_tt = Beta_ttm(t,:)';
    beta_t1t = Beta_tlm(t+1,:)';
    
    P_tt = P_ttm(:,:,t);   
    P_t1t = P_tlm(:,:,t+1); 
    PFP = P_tt*F'/P_t1t; 
    
    % Compute smooth equation (1) and (2)
    beta_tT = beta_tt + PFP*(beta_t1T - beta_t1t);
    P_tT = P_tt + PFP*(P_t1T - P_t1t)*PFP';
    
    % Save
    Beta_tTm(t,:) = beta_tT;
    P_tTm(t,:) = diag(P_tT)';
    
    % Updating
    beta_t1T = beta_tT;
    P_t1T = P_tT;
    t = t - 1;                   
    
end

end
