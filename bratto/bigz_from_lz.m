function [ Z ] = bigz_from_lz( z,z1,z2 )
Z = (z-.5*(z1+z2))/(.5*(z2-z1));
end

