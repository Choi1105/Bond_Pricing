%% makeBeta_ll Fuction
% Unconditional mean
% input: F, Mu (paramters)
% output: beta_ll (Unconditional mean)

function [beta_ll] = makeBeta_ll(F,Mu) 

% Dimension
k = rows(F); 

% Unconditional mean
eyeF = eye(k) - F; 
beta_ll = eyeF\Mu; 

end 