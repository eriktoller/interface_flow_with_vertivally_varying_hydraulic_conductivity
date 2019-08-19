function [ chi ] = chi_of_lz_e( z,nu,z_1,z_2 )
if nu == 0
    chi = (2*z-z_1-z_2)/(z_2-z_1);
else
Z=(2*z-z_1-z_2)/(z_2-z_1);
chi=nu*(Z+sqrt(Z-1)*sqrt(Z+1));
end
end

