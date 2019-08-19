%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface with varying hydraulic conductivity %
% Solution for interface flow                   %
% Erik Toller                                   %
% 2019-05-28                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

% Definition of Varibales
fi0 = 0;
fi1 = 10;
fi2 = 5;
fi3 = 1;
rho_f = 1000;
rho_s = 1025;
fac = 1.2;
lambda = 10;
lambda2 = 10;
lambda3 = 10;
k0 = 1.2*10^-6;
L = 100;
x0 = 0;
x1 = L;

% Constants calculations
alpha = rho_f/(rho_s-rho_f);
kc = 2*k0*lambda^2*((exp(fi1/lambda)-1)+(1/alpha)*(exp(-alpha*fi1/lambda)-1))/((1+alpha)*fi1^2);
kc2 = 2*k0*lambda2^2*((exp(fi2/lambda2)-1)+(1/alpha)*(exp(-alpha*fi2/lambda2)-1))/((1+alpha)*fi2^2);
kc3 = 2*k0*lambda3^2*((exp(fi3/lambda3)-1)+(1/alpha)*(exp(-alpha*fi3/lambda3)-1))/((1+alpha)*fi3^2);

% Uniform flow
Phi0 = Phi_from_fi_dimless(fi0,lambda,alpha)./k0;
Phi1 = Phi_from_fi_dimless(fi1,lambda,alpha)./k0;
Qx0 = (Phi1-Phi0)/L;
Phi02 = Phi_from_fi_dimless(fi0,lambda2,alpha)./k0;
Phi12 = Phi_from_fi_dimless(fi2,lambda2,alpha)./k0;
Qx02 = (Phi12-Phi02)/L;
Phi03 = Phi_from_fi_dimless(fi0,lambda3,alpha)./k0;
Phi13 = Phi_from_fi_dimless(fi3,lambda3,alpha)./k0;
Qx03 = (Phi13-Phi03)/L;

disp('Phi0')
disp(Phi0)


% Tester
disp('Start Check')
inc = 0.001;
[fi_array,Phi_head] = make_fi_array(@(fi)Phi_from_fi_dimless(fi,lambda,alpha)./k0,0,fi1*1.1,inc);
func_fi = @(x)fi_solver(x,@(x)Phi_total(x,Qx0,Phi0),Phi_head, fi_array);
[fi_array2,Phi_head2] = make_fi_array(@(fi)Phi_from_fi_dimless(fi,lambda2,alpha)./k0,0,fi2*1.1,inc);
func_fi2 = @(x)fi_solver(x,@(x)Phi_total(x,Qx02,Phi02),Phi_head2, fi_array2);
[fi_array3,Phi_head3] = make_fi_array(@(fi)Phi_from_fi_dimless(fi,lambda3,alpha)./k0,0,fi3*1.1,inc);
func_fi3 = @(x)fi_solver(x,@(x)Phi_total(x,Qx03,Phi03),Phi_head3, fi_array3);


fi_test(1) = func_fi(x0);
fi_test(2) = func_fi(x1);
fi_check = abs(fi_test/[fi0,fi1]);
if fi_check > 0.99
    disp('Check OK')
else
    disp('Check NOT OK')
end


% Plotting Solver
disp('Start plotting')
xx = linspace(x0-10,x1,100);
fi = func_fi(xx);
fi_2 = func_fi2(xx);
fi_3 = func_fi3(xx);

[h_s] = h_s_from_fi(fi,alpha);
[h_s2] = h_s_from_fi(fi_2,alpha);
[h_s3] = h_s_from_fi(fi_3,alpha);

Phi0_const = Phi_from_fi_const_dimless(fi0,kc,k0,lambda,alpha);
Phi1_const = Phi_from_fi_const_dimless(fi1,kc,k0,lambda,alpha);


fi_const = real(fi_from_Phi_const(xx,@(x)Phi_total(x,Qx0,0),kc./k0,k0,lambda,alpha)) ;
[h_s_const] = h_s_from_fi(fi_const,alpha);
fi_const2 = real(fi_from_Phi_const(xx,@(x)Phi_total(x,Qx02,0),kc2./k0,k0,lambda2,alpha)) ;
[h_s_const2] = h_s_from_fi(fi_const2,alpha);
fi_const3 = real(fi_from_Phi_const(xx,@(x)Phi_total(x,Qx03,0),kc3./k0,k0,lambda3,alpha)) ;
[h_s_const3] = h_s_from_fi(fi_const3,alpha);

% Plotting
figure(1)
hold on
plot(xx./lambda,fi./lambda,'blue -')
plot(xx./lambda,h_s./lambda,'red -')
plot(xx./lambda,fi_const./lambda,'blue -.')
plot(xx./lambda,h_s_const./lambda,'red -.')

plot(xx./lambda2,fi_2./lambda2,'blue -')
plot(xx./lambda2,h_s2./lambda2,'red -')
plot(xx./lambda2,fi_const2./lambda2,'blue -.')
plot(xx./lambda2,h_s_const2./lambda2,'red -.')

plot(xx./lambda3,fi_3./lambda3,'blue -')
plot(xx./lambda3,h_s3./lambda3,'red -')
plot(xx./lambda3,fi_const3./lambda3,'blue -.')
plot(xx./lambda3,h_s_const3./lambda3,'red -.')

title('Interface Flow Comparison, k=k(z) and k=k_{c}')
xlabel('x/\lambda-direction [-]')
ylabel('z/\lambda [-]')

text(xx(70)./lambda,h_s(70)./lambda+1,'k=k(z)/k_0')
text(xx(70)./lambda,h_s_const(70)./lambda-4,'k=k_{c}/k_0')

text(xx(100)./lambda,h_s(100)./lambda,'\phi_1/\lambda = 1')
text(xx(100)./lambda2,h_s2(100)./lambda2,'\phi_1/\lambda = 0.5')
text(xx(100)./lambda3,h_s3(100)./lambda3,'\phi_1/\lambda = 0.1')

axis([min(xx./lambda) max(xx./lambda)+2 min(h_s)./lambda-5 max(fi)./lambda+2])

