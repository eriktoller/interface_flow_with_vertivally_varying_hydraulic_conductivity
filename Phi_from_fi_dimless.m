function [Phi] = Phi_from_fi_dimless(fi,lambda,alpha)
Phi = ( (exp(fi./lambda) - 1) + (1/alpha).*(exp(-alpha.*fi./lambda)-1) );
end

