function [ Omega ] = Omega_lake( chi,a,Q,chi_far,N )
if chi*conj(chi)<0.999
  Omega=complex(NaN,NaN);
else
  if(N==0)
      Omega=0;
  else
      Omega=a(N+1);
      for m=1:N
          mm=N-m+1;
          if mm~=1
            Omega=Omega/chi+a(mm);
          else
            Omega=Omega/chi;
          end
      end
  end
  if Q~=0
    Omega=Omega+Q/(2*pi)*log(chi/abs(chi_far));
  end
end

