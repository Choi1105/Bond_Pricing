function y = lnpdfWishart(Sigma, W, v)
% Compute log pdf of a Wishart distribution.
% Input:
%   Sigma: d x d covariance matrix
%   W: d x d covariance parameter
%   v: degree of freedom
% Output:
%   y: probability density in logrithm scale y=log p(Sigma)
% Written by Mo Chen (sth4nth@gmail.com).
d = length(Sigma);
B = -0.5*v*log(det(W))-0.5*v*d*log(2)-logmvgamma(0.5*v,d);
y = B+0.5*(v-d-1)*log(det(Sigma))-0.5*trace(W\Sigma);

end

function y = logmvgamma(x,d)
% Compute logarithm multivariate Gamma function.
% Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
% log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
% Written by Michael Chen (sth4nth@gmail.com).
s = size(x);
x = reshape(x,1,prod(s));
x = bsxfun(@plus,repmat(x,d,1),(1-(1:d)')/2);
y = d*(d-1)/4*log(pi)+sum(gammaln(x),1);
y = reshape(y,s);
end