function lambda = makeLambda(theta,Spec)

k = Spec.k;
lambda = zeros(k, 1);
lambda(1:k, 1) = theta(Spec.indlambda);

end

