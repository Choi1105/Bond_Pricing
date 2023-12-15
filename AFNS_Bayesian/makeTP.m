% Risk_tsm: T by 1
function [Risk_tsm, Risk_const, Risk_time_varying, Risk_time_varying_Latent, Risk_time_varying_Macro] = makeTP(theta, Xm, B_full, tauj, Spec) 
  
[lambda, ~, G, Omega, ~, GQ] = make_Para(theta, Spec);

[~, Phi_f, Phi_m] = makePhi(G, GQ);
k = rows(lambda);
Omega_ff = Omega(1:k, 1:k);

[Risk_tsm, Risk_const, Risk_time_varying, Risk_time_varying_Latent, Risk_time_varying_Macro] = ...
makeTP_sub_v3_mex(Xm, Omega_ff, B_full, lambda, Phi_f, Phi_m, tauj);   
    
   
end