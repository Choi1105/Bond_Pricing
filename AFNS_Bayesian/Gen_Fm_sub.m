function Fm = Gen_Fm_sub(YCm, ind_basis, a, b)

T = rows(YCm);
k = length(ind_basis);
abar = a(ind_basis);
bbar = b(ind_basis, 1:k);
bbar_inv = inv(bbar);

Fm = (YCm(:, ind_basis) - kron(ones(T,1), abar'))*bbar_inv';

end