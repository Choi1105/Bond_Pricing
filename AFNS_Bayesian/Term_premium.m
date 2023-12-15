function [Risk_tsm_3Dm, Risk_tsm_3D, Risk_const, Risk_time_varying, Risk_time_varying_Latent, Risk_time_varying_Macro] = ...
    Term_premium(post_param, Spec)

Macrom = Spec.Macrom;
mz = cols(Macrom);

n1 = rows(post_param);
tau = Spec.tau; % the number of maturities
ntau = rows(tau);
YCm = Spec.YCm;

T = rows(YCm);

Risk_tsm_3Dm = zeros(n1, T, ntau); % to save risk premium

Risk_const = zeros(T, ntau);
Risk_time_varying_Latent = zeros(T, ntau);
Risk_time_varying_Macro = Risk_time_varying_Latent;

for iter = 1:n1
    
    psi = post_param(iter, :)';
    theta = makeTheta(psi, Spec);
    [Fm, a, b, ~, B_full] = Gen_Fm(theta, Spec);

    if mz > 0
        Xm = [Fm, Macrom];
    else
        Xm = Fm;
    end
%     plot(1200*Fm)
%     Xm = demeanc(Xm);
    %% Term premium computation
    
        for indtau = 2:ntau
            [Risk_tsm_iterj, Risk_const_j, ~, Risk_time_varying_Latent_j, Risk_time_varying_Macro_j] = ...
                makeTP(theta, Xm, B_full, tau(indtau), Spec);
            Risk_tsm_3Dm(iter, :, indtau) = Risk_tsm_iterj;
            Risk_const(:, indtau) = Risk_const(:, indtau) + Risk_const_j;
            Risk_time_varying_Latent(:, indtau) = Risk_time_varying_Latent(:, indtau) + Risk_time_varying_Latent_j;
            Risk_time_varying_Macro(:, indtau) = Risk_time_varying_Macro(:, indtau) + Risk_time_varying_Macro_j;
        end

end

Risk_const = Risk_const/n1;
Risk_time_varying_Latent = Risk_time_varying_Latent/n1;
Risk_time_varying_Macro = Risk_time_varying_Macro/n1;
Risk_time_varying = Risk_time_varying_Latent + Risk_time_varying_Macro;


%% Summary of Output
Risk_tsm_3D = zeros(T, ntau);
for indtau = 1:ntau
    Risk_tsm_3D(:, indtau) = meanc(Risk_tsm_3Dm(:, :, indtau));
end


disp('Term premium computation completes!!!');


end

