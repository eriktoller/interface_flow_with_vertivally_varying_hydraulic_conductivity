function [ z ] = lz_of_chi_e( chi,nu,z_1,z_2)
if nu == 0
    z = .5.*(chi.*(z_2-z_1)+z_2+z_1);
else
Z=.5.*(chi./nu+nu./chi);
z=.5.*(Z.*(z_2-z_1)+z_2+z_1);
end
end

