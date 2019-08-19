function [ Phi ] = Phi_arel( z,zc1,zc2,g0 )
Phi=0;
n = length(zc1);
for j=1:n
    Phi=Phi+g0(j)*Phi_arel_side(z,zc1(j),zc2(j));
end
end

