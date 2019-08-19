function [fi] = fi_from_Phi_const(x,Phi_func,kc,k0,lambda,alpha)
Phi = Phi_func(x);
fi = sqrt(2*Phi*k0*lambda^2/(kc*(1+alpha)));
end

