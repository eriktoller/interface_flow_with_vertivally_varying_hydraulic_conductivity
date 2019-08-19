function [ Phi ] = Phi_arel_side( z,z1,z2)
%Complex potential for the side of an area-element
nu=z2-z1;
Z=(2*z-z1-z2)/nu;
Lsq=nu*conj(nu);
F=-(Z+1)*log(Z+1)+(Z-1)*log(Z-1)+2-2*log(.5*nu);
Phi=Lsq/(64*pi*1i)*(conj(Z)-Z)*(F+conj(F));
end
