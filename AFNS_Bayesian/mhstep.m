% M-H step
function [psi, accept, lnlik0, lnprior0, lnpost0] = mhstep(lnpost, lnlik,lnprior,paramconst,psi,indbj,nu, Spec)

lnlik0 = lnlik(psi, Spec);
lnprior0 = lnprior(psi, Spec);
lnpost0 = lnlik0 + lnprior0;   %  log of denominator
accept = 0;  % = 1 if a proposal is accepted

%% proposal step
[psi1, psimx, Vm, Pm, err1] = Gen_proposal_TaRB(lnpost, paramconst, psi, indbj, nu, Spec);

%% M-H rate step
if err1 < 0.5  % err1 = 0 if Omega is symmetric and positive definite

    % check if psi1 satisfies the parameter constraints
    valid = paramconst(psi1, Spec);
    % if it does not, it is rejected and psi is retained

    if valid > 0.5  % i.e. valid == 1
        % if it does, continue to compute the likelihood and prior at the proposed value */
        lnprior1 = lnprior(psi1, Spec);  %  log prior at psi1
        lnlik1 = lnlik(psi1, Spec);  % log lik at psi1
        lnpost1 = lnlik1 + lnprior1; % log of numerator

        q1 = lnpdfmvt1(psi1(indbj), psimx(indbj), Pm, nu);
        q0 = lnpdfmvt1(psi(indbj), psimx(indbj), Pm, nu);

        loga = lnpost1 + q0 - lnpost0 - q1; % log of probability of acceptance

        accept = log(rand(1,1)) < loga;  % 1 if accepted, and 0 otherwise
        psi = accept*psi1 + (1 - accept)*psi;
        lnpost0 = accept*lnpost1 + (1 - accept)*lnpost0;
        lnlik0 = accept*lnlik1 + (1 - accept)*lnlik0;
        lnprior0 = lnpost0 - lnlik0;

    end
end

end