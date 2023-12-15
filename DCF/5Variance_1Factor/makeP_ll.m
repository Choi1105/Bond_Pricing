%% makeP_ll Fuction
% Unconditional variance
% input: F, Q (paramters)
% output: P_ll (Unconditional variance)

function [P_ll] = makeP_ll(F,Q) 

% Dimension
k = rows(F); 

% Vectorization
Qvec = vec(Q);
kronF = kron(F,F); 
eyekronF = eye(k) - kronF;

% Unconditional variance
P_ll = (eyekronF)\Qvec; 
P_ll = reshape(P_ll,k,k); 
P_ll = (P_ll + P_ll')/2; 

end 