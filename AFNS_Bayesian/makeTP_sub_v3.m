function [Risk_tsm, Risk_const, Risk_time_varying, Risk_time_varying_Latent, Risk_time_varying_Macro] = makeTP_sub_v3(Xm, Omega_ff, B_full, lambda, Phi_f, Phi_m, tauj)

sf = 1200; % scale factor
k = length(Phi_f);
T = rows(Xm);

B_full =  B_full(:, 1:tauj-1); % B_full = 3 by 120

Jensen = 0.5*diag(B_full'*Omega_ff*B_full);

%% 시변하지 않는 부분
Risk_const = sumc(-B_full'*lambda - Jensen)/tauj;
Risk_const = sf*ones(T, 1)*Risk_const; % annual percent로 단위 조정

% disp([meanc(BL*Phi_m) meanc(BL*Phi_f)])
% Phi_f(1,1) = 0;
%% 시변하는 부분
Risk_time_varying_Latent = zeros(T, 1);
Risk_time_varying_Macro = zeros(T, 1);

for t = 1:T
    
    ft = Xm(t, :)';

    RP_TV_Latent = 0;
    RP_TV_Macro = 0;
    for ind_tau = 2:tauj-1
        RP_TV_Latent = RP_TV_Latent - B_full(:, ind_tau)'*Phi_f*ft(1:k);
        RP_TV_Macro = RP_TV_Macro - B_full(:, ind_tau)'*Phi_m*ft(k+1:end);
    end
    
    Risk_time_varying_Latent(t) = RP_TV_Latent;
    Risk_time_varying_Macro(t) = RP_TV_Macro;
    
end

Risk_time_varying_Latent = sf*Risk_time_varying_Latent/tauj;  % annual percent로 단위 조정
Risk_time_varying_Macro = sf*Risk_time_varying_Macro/tauj;

Risk_time_varying = Risk_time_varying_Macro + Risk_time_varying_Latent;
Risk_tsm = Risk_const + Risk_time_varying;

end