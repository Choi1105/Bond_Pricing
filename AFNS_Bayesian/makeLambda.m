function lambda = makeLambda(theta,Spec)

k = Spec.k;
lambda = theta(Spec.ind_lambda);

% Lambda중에 1,2번째만 추정하고 3번째는 0으로 고정하는 경우에 사용 
if length(lambda) < k
  lambda = [lambda; zeros(k-rows(lambda),1)];
end

end

