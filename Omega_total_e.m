function [ Omega ] = Omega_total_e( z,m_not,z_1,z_2,nu,a,Q,chi_far,M,m_coef,C,zc1,zc2,g0)
Omega=Phi_arel( z,zc1,zc2,g0 )+C;

if M>0
    for mm=1:M
        if mm ~= m_not
            chi=chi_of_lz_e(z,nu(mm),z_1(mm),z_2(mm));
            Omega=Omega+Omega_lake(chi,a(mm,:),Q(mm),chi_far(mm),m_coef(mm));
        end
    end
end
end

