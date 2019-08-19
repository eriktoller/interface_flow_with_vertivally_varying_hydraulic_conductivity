function [a,Q,C] = solve_lakes_e(Phi0,Phi_lake,M,N,m,z_1,z_2,nu,chi_far,zc1,zc2,g0,z_ref)
error=1e6;
NIT=0;
C = 0;
m_max = max(m);
a=zeros(M,m_max+1);
Q_old=zeros(1,M);
erQmax=0;
eramax=0;
eramaxr=0;
a_old=a;
Qsum=0;
Q=Q_old;
asum=zeros(1,M);
NITmax = 2000;
N2 = round(N/10);
AM=get_AMQ_e(M,N2,z_1,z_2,nu,chi_far,z_ref);

while error>1e-6 && NIT<NITmax
        Q=solve_Q_e( Phi_lake,M,N2,m,a,z_1,z_2,nu,chi_far,AM,zc1,zc2,g0,Phi0,z_ref);
        C=Q(M+1);
    for mm=1:M
        a(mm,:)=-conj(Cauchy_integral(N2,m(mm),m_max,@(chi)lz_of_chi_e(chi,nu(mm),z_1(mm),z_2(mm))...
            ,@(z)Omega_total_e( z,mm,z_1,z_2,nu,a,Q,chi_far,M,m,C,zc1,zc2,g0)));
        erQ=abs(Q(mm)-Q_old(mm));
        Qsum=Qsum+abs(Q(mm));
        if erQ>erQmax
            erQmax=erQ;
        end
    end
    for mm=1:M
        for kk=1:m+1
            era=abs(a(mm,kk)-a_old(mm,kk));
            asum(mm)=asum(mm)+abs(a(mm,kk));
            if era>eramax
                eramax=era;
            end
        end
        erar=M*eramax/asum(mm);
        if erar>eramaxr
            eramaxr=erar;
        end
    end
    disp('Discharge of lake & change')
    disp([Q(26),Q(26)-Q_old(26)])
    NIT=NIT+1;
    disp('NIT')
    disp(NIT)
    erQmaxr=M*abs(erQmax/Qsum);
    if erQmaxr>eramaxr
        error=erQmaxr;
    else
        error=eramaxr;
    end
    disp('error Q & error a')
    disp([erQmaxr,eramaxr])
    a_old=a;
    Q_old=Q;
    Qsum=0;
    asum=zeros(1,M);
    erQmax=0;
    eramax=0;
    eramaxr=0;
end
end

