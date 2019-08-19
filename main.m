%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface with varying hydraulic conductivity %
% Solution for Bratto island                    %
% Erik Toller                                   %
% 2019-05-28                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
format long
tic
load geometry
%% Definition of Varibales
% Aquifer
fi0 = 0;
rho_f = 1000;
rho_s = 1025;
lambda = 140;
k0 = 0.0047;
eta = 200;
H_s = fi0;

% Hydraulic Conductivty Calculations
alpha = rho_f/(rho_s-rho_f);

% Reference point
z_ref = complex(661500,6423500);
Phi0 = Phi_from_fi(fi0,k0,lambda,alpha);

% Wells
rw = bratto_wells(:,3)./1000;
fi_bottom_well = bratto_wells(:,4);
fi_max_well = fi_bottom_well./alpha;
zw = complex(bratto_wells(:,1),bratto_wells(:,2));
nw = length(zw);


% Analytical Elements (island boundary)
z_1 = complex(bratto_boundary(:,1),bratto_boundary(:,2));
z_2 = z_1;
z_2(1) = [];
z_2(end+1) = z_1(1);
z_1 = [z_1;complex(bratto_lake(:,1),bratto_lake(:,2))];
z_2 = [z_2;complex(bratto_lake(:,3),bratto_lake(:,4))];

% adding the wells (as analytical elements)
num = length(z_1);
for ii = 1:nw
    z_1(num+ii) = zw(ii)+complex(rw(ii),0);
    z_2(num+ii) = zw(ii)-complex(rw(ii),0);
end

% Pre-calculations and defentions for iterative sovler
M = length(z_1);
nu = [linspace(1,1,num-2),0.3,1,linspace(0,0,nw)];
fi_lake = [linspace(H_s,H_s,num-2),60,10,fi_max_well'];
m = 40;
m_coef = [linspace(m,m,num),linspace(1,1,nw)];
N = 800;

chi_far = zeros(1,M);
for ii = 1:M
    chi_far(ii) = chi_of_lz_e(z_ref,nu(ii),z_1(ii),z_2(ii));
end

Phi_lake = zeros(1,M); 
for ii = 1:M
    Phi_lake(ii) = Phi_from_fi(fi_lake(ii),k0,lambda,alpha);
end

% Infiltration
range = 1:length(bratto_infiltration);
zc1 = complex(bratto_infiltration(range,1),bratto_infiltration(range,2));
zc2 = complex(bratto_infiltration(range,3),bratto_infiltration(range,4));
g0 = bratto_infiltration(:,5)*20;

%% Coefficient solver
[a,Q,C] = solve_lakes_e(Phi0,Phi_lake,M,N,m_coef,z_1,z_2,nu,chi_far,zc1,zc2,g0,z_ref);

%% Head functions
inc = 0.001;
[fi_array,Phi_head] = make_fi_array(@(fi)Phi_from_fi(fi,k0,lambda,alpha),min(fi_lake),eta,inc);

func_fi = @(z)fi_solver(z,@(z)real(Omega_total_e( z,0,z_1,z_2,nu,a,Q,chi_far,M,m_coef,C,zc1,zc2,g0)),...
    Phi_head, fi_array);

func_Phi = @(z)real(Omega_total_e( z,0,z_1,z_2,nu,a,Q,chi_far,M,m_coef,C,zc1,zc2,g0));

func_Omega = @(z) Omega_total_e( z,0,z_1,z_2,nu,a,Q,chi_far,M,m_coef,C,zc1,zc2,g0);


%% Plotting data
xfrom = min(real(z_1))-100;
xto = max(real(z_1))+100;
Nx = 100;
yfrom = min(imag(z_1))-100;
yto = max(imag(z_1))+100;
Ny = Nx;
nint = 30;

%% Plotting
%% Discharge Potential

disp('Plotting flow net')
figure;
ContourMe_flow_net(xfrom,xto,Nx,yfrom,yto,Ny,func_Phi,nint);

legend off
axis equal
hold on

ellipse_c(z_1,z_2,nu,'black');
ellipse_c(zc1,zc2,linspace(1,1,length(zc1)),'red');
plot(z_ref,'blue d')
plot(zw,'blue *')
title('Discharge Potential')
xlabel('E-coordinate [-]')
ylabel('N-coordinate [-]')
disp('done')
  
%% Head contour
disp('Plotting head contour')
figure;
fi_grid = ContourMe_fi(xfrom,xto,Nx,yfrom,yto,Ny,func_fi,nint);

legend off
axis equal
hold on

ellipse_c(z_1,z_2,nu,'black');
ellipse_c(zc1,zc2,linspace(1,1,length(zc1)),'red');
plot(z_ref,'blue d')
plot(zw,'blue *')
title('Hydraulic Head')
xlabel('E-coordinate [-]')
ylabel('N-coordinate [-]')
disp('done')

%% Interface contour
disp('Plotting interface contour')
figure;
h_s_grid = h_s_from_fi(fi_grid,alpha);
X = linspace(xfrom, xto, Nx);
Y = linspace(yfrom, yto, Ny);
contour(X, Y,h_s_grid,nint);
colorbar

legend off
axis equal
hold on

ellipse_c(z_1,z_2,nu,'black');
ellipse_c(zc1,zc2,linspace(1,1,length(zc1)),'red');
plot(z_ref,'blue d')
plot(zw,'blue *')
title('Interface Head')
xlabel('E-coordinate [-]')
ylabel('N-coordinate [-]')
disp('done')

%% Hydraulic conductivity plot
disp('Plotting hydraulic conductivity')
fi_plot = linspace(-100,100,N);
k_plot = k_from_fi(fi_plot,k0,lambda);
figure;
hold on
plot(k_plot,fi_plot,'blue')
grid minor
title('Hydraulic Conductivity')
xlabel('Hydraulic conductivity [m/day]')
ylabel('Depth [m]')
disp('done')


%% Discharges (wells)

figure;
hold on
grid minor

bar(Q(28:60))
title('Maximum Discharge of Wells')
ylabel('Discharge [m^3/day]')
xlabel('Well number')
%% Map

figure;
axis equal
hold on

ellipse_c(z_1,z_2,nu,'black');
ellipse_c(zc1,zc2,linspace(1,1,length(zc1)),'red');

xlabel('E-coordinate [-]')
ylabel('N-coordinate [-]')
title('Map Brattö')

plot(zw,'* blue')
formatSpec = ' %d ';
for ii = 1:nw
str = sprintf(formatSpec,ii);
tx = text(real(zw(ii)),imag(zw(ii)),str);
if mod(ii,2)
    tx.HorizontalAlignment = 'right';
else
    tx.HorizontalAlignment = 'left';
end
tx.FontSize = 8;
end

tx = text(real(z_1(11))-480,imag(z_2(11)),'Island boundary \phi = 0 m \rightarrow');
tx.FontSize = 9;
tx = text(real(z_1(27))-300,imag(z_2(27))+200,'Ditch \phi = 10 m \rightarrow');
tx.FontSize = 9;
tx = text(real(z_1(26))-400,imag(z_2(26)),'Lake \phi = 60 m \rightarrow');
tx.FontSize = 9;
tx = text(real(zc1(2)),imag(zc2(2)),'\gamma_1');
tx.FontSize = 10;
tx.Color = 'red';
tx = text(real(zc1(10)),imag(zc2(10))+20,'\gamma_2');
tx.FontSize = 10;
tx.Color = 'red';
tx = text(real(zc1(19)),imag(zc2(19)),'\gamma_3');
tx.FontSize = 10;
tx.Color = 'red';
tx = text(real(zc1(21)),imag(zc2(21)),'\gamma_4');
tx.FontSize = 10;
tx.Color = 'red';

%% Check
fi_ref = func_fi(z_ref);
fi_boundary = zeros(M,1);
for ii = 1:M
    fi_boundary(ii) = func_fi(lz_of_chi_e_array(1,nu(ii),z_1(ii),z_2(ii)));
end

disp('Head at reference point')
disp(fi_ref)
disp('Head at boundaries')
disp(fi_boundary)

disp('---------------')
toc
disp('---------------')