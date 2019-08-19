function [ AM ] = get_AMQ_e(M,N,z_1,z_2,nu,chi_far,z_ref )
AM=zeros(M+1,M+1);
for jj=1:M
    for mm=1:M
        a=real(Cauchy_integral(N,0,0,@(z)lz_of_chi_e(z,nu(jj),z_1(jj),z_2(jj)),...
            @(z)log(chi_of_lz_e(z,nu(mm),z_1(mm),z_2(mm))/abs(chi_far(mm)))/(2*pi)));
        AM(jj,mm)=a(1);
    end
    AM(jj,M+1)=1;
end

for jj=1:1
    for mm=1:M
        a=real(log(chi_of_lz_e(z_ref,nu(mm),z_1(mm),z_2(mm))/abs(chi_far(mm)))/(2*pi));
        AM(M+jj,mm)=a(1);
    end
    AM(M+jj,M+1)=1;
end
end

