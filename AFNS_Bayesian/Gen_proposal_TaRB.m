function [psi1, psimx, V, Vinv, err] = Gen_proposal_TaRB(lnpost, paramconst, psi,indbj,nu, Spec)

% psi0 = make_Omegahat_MH(psi, Spec);
psi0 = psi;
maxiter = 2; % maximum iteration in Newton step
IT = 0.00001;
trf = 0.0001;   % temperature reduction factor
SL = 0;     % stage length increment
IM = 20;     %length of first stage
cs = 50;
NS = 1;    % number of SA stages
eps = 1e-6;      % desired precision of maximum
printi = 0;         % 1 to print interim output, 0 otherwise
mr = 400;    % number of rejection
co = 1;   % 1 if constrained optimization, and 0 if unconstrained optimization
SF = 100*ones(rows(indbj),1);    % scale factor

[psimx, ~, ~, V, Vinv] = SA_CRK2(lnpost,paramconst,SF,psi0,indbj,trf,IT,SL,cs,IM,NS,mr,eps,maxiter,co,printi,Spec);

V = diag(diag(V));

% proposal from t-distribution
psi1 = psi;
[t_val, err] = randmvt(nu, psimx(indbj), V);
psi1(indbj) = t_val; % proposal


end