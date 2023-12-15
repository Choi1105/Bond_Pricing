%% Compute log lnlikelihood
% input: psi, Sn (paramters, structure parameters)
% output: lnL (log likelihood)

function [lnL] = lnlik(psi,Sn)

% Data
ym = Sn.data; 

% Parameter transform / state space parameter form
theta = maketheta(psi,Sn);
[C,H,R,Mu,F,Q] = makePara(theta,Sn);

% Compute likelihood by Kalman Filter
[lnL, ~, ~] = KM_filter(C,H,R,Mu,F,Q,ym);

end