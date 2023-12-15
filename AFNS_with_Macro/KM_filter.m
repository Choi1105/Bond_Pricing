% Kalman filter

function [Beta_ttm, P_ttm, Fittedm, Residm] = KM_filter(a, b, G, Q, Sigma, ym)

Null_a = zeros(2,1);
a_New = [a;Null_a];

Null_M = zeros(10,2);
Null_b = zeros(2,5);
b_ew = [b Null_M];
b_New = [b, zeros(10,2) ; zeros(2,3), diag(ones(2,1))];

Sigma(11:12,11:12) = zeros;

[T,ntaum] = size(ym);
k = rows(G);
zerom = zeros(ntaum,1);
f_ll = zeros(k,1);
P_ll = makeR0(G,Q);
ek = eye(k);
lnL = 0;

Beta_ttm = zeros(T,k);
P_ttm = zeros(T,k);
Fittedm = zeros(T,ntaum);
Residm = zeros(T,ntaum);

for t = 1:T

    f_tl = G*f_ll;
    P_tl = G*P_ll*G' + Q;
    var_tl = Sigma + b_New*P_tl*b_New';
    var_tl = 0.5*(var_tl + var_tl');

    e_tl = ym(t,:)' - a_New - b_New*f_tl;
    Kalgain = P_tl*b_New';
    Kalgain = Kalgain/var_tl;
    f_tt = f_tl + Kalgain*e_tl;
    P_tt = ek - Kalgain*b_New;
    P_tt = P_tt*P_tl;

    lnL = lnL + lnpdfmvn(e_tl,zerom,var_tl);

    f_ll = f_tt;
    P_ll = P_tt;
    
    Beta_ttm(t,:) = f_tt';
    P_ttm(t,:) = diag(P_tt)';
    Pred_Ym = a_New + b_New*f_tt;
    Fittedm(t,:) = Pred_Ym';
    Residm(t,:) = ym(t,:)' - Pred_Ym;
end

end
