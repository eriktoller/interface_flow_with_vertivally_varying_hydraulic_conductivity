function [k] = k_from_fi(fi,k0,lambda)
k = k0.*exp(fi./lambda);
end

