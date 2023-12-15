% Parameters in S-S form
function [Lam,Sigma,mu,G,L,Gamma] = makePara(theta, Sn) 

tau = Sn.tau;
lambda = Sn.lambda;

% M.E.
ntau = rows(Sn.tau);
Lam = ones(ntau, 3);
Lam(:,2) = (ones(ntau,1)-exp(-tau*lambda))./(tau*lambda);
Lam(:,3) = ((ones(ntau,1)-exp(-tau*lambda))./(tau*lambda)) - exp(-tau*lambda);

Sigma = diag(theta(Sn.indSig));

% T.E.
mu = theta(Sn.indMu); 
G = reshape(theta(Sn.indG),Sn.k,Sn.k);

% Omega = L * Gamma * L' (Spectural Decomposition)

% L     = [  Sig(1)     0         0     ;...
%              0      Sig(2)      0     ;...
%              0        0       Sig(3)  ]

% Gamma     = [  1       gamma(1)   gamma(2) ;...
%              gamma(1)      1      gamma(3) ;...
%              gamma(2)    gamma(3)     1    ]

sig2 = theta(Sn.indLmatrix);
L = [sig2(1)    0        0     ;...
       0      sig2(2)    0     ;...
       0        0      sig2(3) ];

Rho = theta(Sn.indGamma);
Gamma = [  1     Rho(1)  Rho(2) ;...
         Rho(1)    1     Rho(3) ;...
         Rho(2)  Rho(3)    1    ];

end