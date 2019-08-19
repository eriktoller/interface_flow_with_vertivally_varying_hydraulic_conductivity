function [k] = k_from_z(z,k0,lambda)
k = k0.*exp(z./lambda);
end

