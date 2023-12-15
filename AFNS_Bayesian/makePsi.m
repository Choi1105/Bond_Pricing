function psi = makePsi(theta, Spec)

% scale factor를 곱해준 parameter psi
% 사전정보는 psi에 대해서 설정해줌
% theta는 scale이 너무 작아서 prior를 정해주는것이 어렵다.

psi = theta;
sf_Omega = Spec.sf_Omega;
sf_Sig = Spec.sf_Sig;
sf_lambda = Spec.sf_lambda;
psi(Spec.ind_V) = sf_Omega*theta(Spec.ind_V);
psi(Spec.ind_Omega_MZ) = (sf_Omega^2)*theta(Spec.ind_Omega_MZ);
psi(Spec.ind_lambda) = sf_lambda*theta(Spec.ind_lambda);
psi(Spec.ind_Sig) = (sf_Sig^2)*theta(Spec.ind_Sig);
 
end