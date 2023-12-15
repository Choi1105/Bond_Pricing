function [abar, bbar, Abar, Bbar] = makeABbar(delta, G, L, Omega, beta, GQ, lambda, tau, k)
%% Notice %%
% This function generates A(tau), B(tau), a(tau), and b(tau).

% We do not utilize the parameters G and Omega in this function.

%% Input parameter %% 
% delta = Mean of the Short Rate (Dimension : 1 by 1 scalar)

% L   =  Cholesky Component of Omega (Dimension : k by k Lower Triangular Matrix)

% beta =  Restriction matrix (1 1 0)' (Dimension : 3 by 1 vector)

% GQ = Measure the Persistent B(tau) (Dimension : k by k)

% lambda = Risk Aversion Parameter (Dimension : k by 1)

%% Output parameter %% 
% Abar = delta + Abar(tau-1) - 0.5*Bbar(tau-1)'*L'*L*Bbar(tau-1) - Bbar(tau-1)'*L*lambda (Dimension : 1 by 1 scalar)

% Bbar = beta + GQ'Bbar(tau-1)  (Dimension : 3 by 1 vector)

% abar = Abar(tau) / tau (Dimension : 1 by 1 )

% bbar = Bbar(tau) / tau (Dimension : 3 by 1 )

%% Appendix %%
% Our final goal is to compute y(tau, t) using the formula y(tau,t) = abar + bbar' * f(t)

% However, in this specific algorithm, particularly within the KM_filter function, we caculate y(tau, t) as abar + bbar * f(t)

% Therefore, in the final step, save the transpose of bbar as bbar =  bbar'

%% Algorithm Code %%
% Number of tau
N = maxc(tau);

% Pre-allocation
Abar = zeros(N,1);
Bbar = zeros(k,N);

abar = zeros(length(tau),1);
bbar = zeros(k,length(tau));

% Initial Parameter value 
Abar(1,1) = delta;      % Dimension : 1 by 1
Bbar(:, 1) = beta;       % Dimension : 3 by 1

% Calculate the Bbar (Exogenous Variable)
for i = 2:N
    Bbar(:, i) = beta + GQ'*Bbar(:,i-1)
end

% Calculate the bbar
for v = 1:length(tau)
    bbar(:,v) = Bbar(:,tau(v))/tau(v);
end

% Calculate the Atau
for z = 2:N
    Abar(z) = delta + Abar(z-1) - 0.5*Bbar(:,z-1)'*L'*L*Bbar(:,z-1) - Bbar(:,z-1)'*L*lambda;
end

% Calculate the atau and btau
for v = 1:length(tau)
    abar(v) = Abar(tau(v))/tau(v);
end

% Final Step (Transpose the bbar)
bbar = bbar';

end

