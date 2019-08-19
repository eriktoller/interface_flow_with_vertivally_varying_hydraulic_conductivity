function [ z ] = lz_of_chi_e_array( chi,nu,z_1,z_2)
z = zeros(1,length(nu));
for ii = 1:length(nu)
    nu0 = nu(ii);
    z_10 = z_1(ii);
    z_20 = z_2(ii);
if nu0 == 0
    z(ii) = .5.*(chi.*(z_20-z_10)+z_20+z_10);
else
Z=.5.*(chi./nu0+nu0./chi);
z(ii)=.5.*(Z.*(z_20-z_10)+z_20+z_10);
end
end
end

