function [ chi ] = chi_of_lz_ls( z,z_1,z_2 )
Z=(2*z-z_1-z_2)/(z_2-z_1);
chi=Z+sqrt(Z-1)*sqrt(Z+1);
end

