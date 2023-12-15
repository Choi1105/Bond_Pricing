function [Gam_til] = makeGam_til(Gam)

 Gam_til = cholmod(Gam)';

end