function lnL = Kalman(a, b, G, Q, Sigma, ym)
Null_a = zeros(2,1);
a_New = [a;Null_a];

Null_M = zeros(10,2);
Null_b = zeros(2,5);
b_ew = [b Null_M];
b_New = [b, zeros(10,2) ; zeros(2,3), diag(ones(2,1))];

Sigma(11:12,11:12) = zeros;

[T,ntaum] = size(ym);
km = rows(G);
zerom = zeros(ntaum,1);
f_ll = zeros(km,1);
P_ll = makeR0(G,Q);
ek = eye(km);
lnL = 0;

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
      
 end



end


