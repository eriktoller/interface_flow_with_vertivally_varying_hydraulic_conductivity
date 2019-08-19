function [Phi] = Phi_from_fi(fi,k0,lambda,alpha)
Phi = k0.*lambda.^2.*( (exp(fi./lambda) - 1) + (1./alpha).*(exp(-alpha.*fi./lambda)-1) );
end

