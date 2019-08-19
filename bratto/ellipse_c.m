function [] = ellipse_c( z_1,z_2,nu,plot_spec )
for ii = 1:length(nu)
    z1_temp = z_1(ii);
    z2_temp = z_2(ii);
    nu_temp = nu(ii);
    theta=0;
    N=100;
    x=zeros(1,N);
    y=zeros(1,N);
    deltheta=2*pi/N;
    for i=1:N
        theta=theta+deltheta;
        chi=exp(1i*theta);
        z=lz_of_chi_e(chi,nu_temp,z1_temp,z2_temp);
        x(i)=real(z);
        y(i)=imag(z);
    end
    plot(x,y,plot_spec );
end
end

