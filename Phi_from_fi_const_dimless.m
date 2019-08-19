function [Phi] = Phi_from_fi_const_dimless(fi,kc,k0,lambda,alpha)
Phi = .5.*(kc./k0).*(fi./lambda).^2.*(1+alpha);
end

