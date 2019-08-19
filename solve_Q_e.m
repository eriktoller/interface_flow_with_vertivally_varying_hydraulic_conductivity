function [ Q ] = solve_Q_e( Phi_lake,M,N,m_coef,a_in,z_1,z_2,nu,chi_far,AM,zc1,zc2,g0,Phi0,z_ref )
KN=zeros(M+1,1);
Q0=zeros(1,M);
C = 0;
for jj=1:M
    a=real(Cauchy_integral(N,0,0,@(chi)lz_of_chi_e(chi,nu(jj),z_1(jj),z_2(jj)),...
        @(z)Omega_total_e( z,0,z_1,z_2,nu,a_in,Q0,chi_far,M,m_coef,C,zc1,zc2,g0)));
    KN(jj,1)=Phi_lake(jj)-a(1);
end

KN(M+1,1)=Phi0-real(Omega_total_e( z_ref,0,z_1,z_2,nu,a_in,Q0,chi_far,M,m_coef,C,zc1,zc2,g0));

Q=AM\KN;


end

