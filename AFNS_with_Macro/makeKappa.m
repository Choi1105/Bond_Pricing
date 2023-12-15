function GQ = makeKappa(kappa)

GQ = zeros(3, 3);
GQ(1,1) = 1;
GQ(2,2) = exp(-kappa);
GQ(2,3) = kappa*GQ(2,2);
GQ(3,3) = GQ(2,2);

end
