function [ z ] = lz_of_chi_ls( chi,z_1,z_2)
Z=.5*(chi+1/chi);
z=.5*(Z*(z_2-z_1)+z_2+z_1);
end

