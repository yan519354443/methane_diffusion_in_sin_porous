%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%
clear
ny = 100;
nx = 1.8*ny;
Re = 60;
U = .1;
%%%%%%%%%%%%%%%% IMMERSED OBJECT %%%%%%%%%%%%%%%
%SINUISOIDAL
amplitude=.3*ny; %SIN AMPLITUDE
k=2;             %# of WAVES
[y,x] = meshgrid(1:ny,1:nx); 
obst=zeros(nx,ny);
obst=(y>ny/2+ny/2+amplitude*sin(2*pi*(x-1)/((nx-1)/k))) | (y<=ny/2-ny/2-amplitude*sin(2*pi*(x-1)/((nx-1)/k)));
obst(:,[1,ny])=1;
immersed = find(obst);
figure;imshow(obst)

%%%%%%%%%%%%%%%% real diameter of path %%%%%%%%%%%%%%%
width_real_each_grid = 10*10^(-9)/(ny-2);  %100纳米下，每个格子的宽度
pore_width = (sum(obst ~= 1, 2)-2)*width_real_each_grid;

%%%%%%%%%%%%%%%%%%% Molecule and Gas Properties %%%%%%%%%%%%%%%%%%%%%%% 
mass_molecule = 2.658 * 10 ^ (-26);
d = 0.38 * 10 ^ (-9); %diameter of methane molecule
%%%%%%%%%%%%%%%% environment properties %%%%%%%%%%%%%%%
T = 298; % in K
Tc = 190.4; % in K
Gc = -0.4; % the interaction strength of fluid particles
Pc = 4.595;   %kpa
Vm_cr = 98.66;
Zcr = 0.287;
NA = 6.02*10^23;
kb = 1.3e-23; %Boltzmann constant
molecular_weight = 16.04; % g/mol
%%%%%%%%%%%%%%%% OUTER BOUNDS %%%%%%%%%%%%%%%%%
L = ny-2; 
H = y-nx/ny;
u = 4*U/(L^2)*(H*L-H.^2);
v = zeros(nx,ny);
%%%%%%%%%%%%%%%% INSTANTIATION %%%%%%%%%%%%%%%%%
%tau = 0.5+3*U*ny/Re;
w = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
oppbb=[1 4 5 2 3 8 9 6 7];
oppmb=[1 4 5 2 3 9 8 7 6]; %镜面反弹
c=[cx;cy]'; 
% u = zeros(ny,nx);
% v = zeros(ny,nx);
%rho = ones(ny,nx);
rho = 1;
feq=zeros(9,nx,ny);
for i=1:9 
    feq(i,:,:)=w(i)*rho.*(1+3*(c(i,1)*u+c(i,2)*v)+ 9/2*(c(i,1)*u+c(i,2)*v).^2 - 3/2*(u.^2+v.^2)); 
end
f = feq;


iter=0;



fnx = zeros(size(f));
%%%%%%%%%%%%%%%%%%%%%%%% START ITERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(iter<30)
uprev=reshape(u,nx,ny);
vprev=reshape(v,nx,ny);
%%%%%%%%%%%%%%%%%%%%%%%% MACRO VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
rho = sum(f);
% rhou=f(1,:,:)*c(1,1);
% for i=2:9
%     rhou=rhou+f(i,:,:)*c(i,1);
% end;
% rhov=f(1,:,:)*c(1,2);
% for i=2:9
%     rhov=rhov+f(i,:,:)*c(i,2);
% end;
% u=squeeze(rhou./rho);
% v=squeeze(rhov./rho);
% rho=squeeze(rho);    
u=reshape((cx * reshape(f,9,nx*ny)),1,nx,ny)./rho;

v=reshape((cy * reshape(f,9,nx*ny)),1,nx,ny)./rho;

%%%%%%%%%%%%%%%%%%% Effective Knusen Number,  Phi and Pressure %%%%%%%%%%%%%%%%%%%%%%% 
kne = eff_Kn_cal(rho, mass_molecule, d, pore_width);
phi = phi_cal(rho, T, Tc,Gc);
pressure = pressure_cal(rho, phi,Gc);

%%%%%%%%%%%%%%%%%%% number density and three mechanisms of diffusion %%%%%%%%%%%%%%%%%%%%%%% 
[R, a , b] = eos_parameters_cal(Tc, Pc, Vm_cr, Zcr);
V = molar_volume_cal(pressure, R, a, b, T);  %v in m^3/kg
number_density = NA/V;
md = molecular_diffution(number_density,kb,mass_molecule,T, d);
kd = knusen_diffusion(rho, R, T, molecular_weight, pore_width);
td = transition_diffusion(md, kd);
%%%%%%%%%%%%%%%%%%% tau_alpha %%%%%%%%%%%%%%%%%%%%%%% 
delta_t = 10e-9;
tau_alpha = tau_alpha_cal(kne,md,kd,td,delta_t);

%%%%%% BOUNDARIES


% INLET zou-he boundary 

%{
H =[2:(ny-1)]-1.5;

u(:,1,2:(ny-1))=4*U/(L^2)*(H*L-H.^2)*10e-9;

v(:,1,[2:(ny-1)])=0;
rho(:,1,[2:(ny-1)])=1./(1-u(:,1,[2:(ny-1)])).*(sum(f([1,3,5],1,[2:(ny-1)]))+2*sum(f([4,7,8],1,[2:(ny-1)])));
%}
rho(:,1,[2:(ny-1)])=2;
u(:,1,[2:(ny-1)])=1 - 1./(rho(:,1,[2:(ny-1)])).*(sum(f([1,3,5],1,[2:(ny-1)]))+2*sum(f([4,7,8],1,[2:(ny-1)])));

v(:,1,[2:(ny-1)])=0;
%rho(:,1,[2:(ny-1)])=50/(iter+1);
f(2,1,[2:(ny-1)])=f(4,1,[2:(ny-1)])+2/3*rho(:,1,[2:(ny-1)]).*u(:,1,[2:(ny-1)]); 


f(6,1,[2:(ny-1)])=f(8,1,[2:(ny-1)])+1/2*(f(5,1,[2:(ny-1)])-f(3,1,[2:(ny-1)]))+ 1/2*rho(:,1,[2:(ny-1)]).*v(:,1,[2:(ny-1)]) ...
                                + 1/6*rho(:,1,[2:(ny-1)]).*u(:,1,[2:(ny-1)]); 
%disp(f(6,1,[2:(ny-1)]))
f(9,1,[2:(ny-1)])=f(7,1,[2:(ny-1)])+1/2*(f(3,1,[2:(ny-1)])-f(5,1,[2:(ny-1)]))-1/2*rho(:,1,[2:(ny-1)]).*v(:,1,[2:(ny-1)]) ...
                                + 1/6*rho(:,1,[2:(ny-1)]).*u(:,1,[2:(ny-1)]); 

% OUTLET
u(:,nx,[2:(ny-1)])=-1+1./(rho(:,nx,[2:(ny-1)])).*(sum(f([1,3,5],nx,[2:(ny-1)]))+2*sum(f([2,6,9],nx,[2:(ny-1)])));

v(:,nx,[2:(ny-1)])=0;
f(4,nx,[2:(ny-1)])=f(2,nx,[2:(ny-1)])-2/3*rho(:,nx,[2:(ny-1)]).*u(:,nx,[2:(ny-1)]); 

f(8,nx,[2:(ny-1)])=f(6,nx,[2:(ny-1)])+1/2*(f(3,nx,[2:(ny-1)])-f(5,nx,[2:(ny-1)]))- 1/2*rho(:,nx,[2:(ny-1)]).*v(:,nx,[2:(ny-1)]) ...
                                  - 1/6*rho(:,nx,[2:(ny-1)]).*u(:,nx,[2:(ny-1)]); 

f(7,nx,[2:(ny-1)]) = f(9,nx,[2:(ny-1)]) + 1/2*(f(5,nx,[2:(ny-1)])-f(3,nx,[2:(ny-1)]))+ 1/2*rho(:,nx,[2:(ny-1)]).*v(:,nx,[2:(ny-1)]) ...
                                  - 1/6*rho(:,nx,[2:(ny-1)]).*u(:,nx,[2:(ny-1)]); 

%{
for i = 1:9
    f(i,1,[2:(ny-1)]) = w(i)*(rho(:,1,[2:(ny-1)])+3*(u(:,1,[2:(ny-1)])+v(:,1,[2:(ny-1)])))+f(i,2,[2:(ny-1)]) - feq(i,2,[2:(ny-1)]);
    f(i,nx,[2:(ny-1)]) = w(i)*(rho(:,nx,[2:(ny-1)])+3*(u(:,nx,[2:(ny-1)])+v(:,nx,[2:(ny-1)])))+f(i,nx-1,[2:(ny-1)]) - feq(i,nx-1,[2:(ny-1)]);
end
%}
%{
%%%  diffusion boundary %%% refer to https://github.com/temmy222/Lattice-Boltzmann-Method-with-Python/blob/master/d2q9%20diffusion.py
%inlet
f(2,1,[2:(ny-1)])=1/9*rho(:,1,[2:(ny-1)]) + 1/9*rho(:,1,[2:(ny-1)]) - f(4,1,[2:(ny-1)]);
f(6,1,[2:(ny-1)])=1/36*rho(:,1,[2:(ny-1)]) + 1/36*rho(:,1,[2:(ny-1)]) - f(8,1,[2:(ny-1)]);
f(9,1,[2:(ny-1)])=1/36*rho(:,1,[2:(ny-1)]) + 1/36*rho(:,1,[2:(ny-1)]) - f(7,1,[2:(ny-1)]);
%outlet
f(5,nx,[2:(ny-1)]) = - f(3,nx,[2:(ny-1)]);
f(8,nx,[2:(ny-1)]) = - f(6,nx,[2:(ny-1)]);
f(9,nx,[2:(ny-1)]) = - f(7,nx,[2:(ny-1)]);
f(4,nx,[2:(ny-1)]) = - f(2,nx,[2:(ny-1)]);
f(1,nx,[2:(ny-1)]) = 0.0;
%}



%%%%%%%%%%%%%%%%%%%% F EQIILIBRIUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%collOperator=(1-1/tau)*f+1/tau*feq; 
for i=1:9
feq(i,:,:)=w(i)*rho.*(1+3*(c(i,1)*u+c(i,2)*v)+ 9/2*(c(i,1)*u+c(i,2)*v).^2 - 3/2*(u.^2+v.^2)); 

    for rows = 1:nx
        for cols = 1:ny
            fnx(i,rows,cols)=f(i,rows,cols)-(1/(tau_alpha(1,rows,cols)))*(f(i,rows,cols)-feq(i,rows,cols));
        end
    end
end
%%%%%%%%%%%%%%%%%%%% COLLISION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOUNCE-BACK

fnx_oppbb = fnx(:,:,:);
fnx_oppmb = fnx(:,:,:);
for i=1:9
fnx_oppbb(i,immersed)=f(oppbb(i),immersed);
fnx_oppmb(i,immersed)=f(oppmb(i),immersed);
end

r_coefficient_value = r_coefficient(kne, pore_width,width_real_each_grid);
[p,m,n] = size(kne);

for rows = 1:m
    for cols = 1:n
        fnx(:, rows, cols) = r_coefficient_value(1,rows, cols)*fnx_oppbb(:,rows, cols) + (1 - r_coefficient_value(1,rows, cols))...
                *fnx_oppmb(:,rows, cols);
    end
end

%%%%%%%%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:9
f(i,:,:)=circshift(fnx(i,:,:),[0,cx(i),cy(i)]);
end
%%%%%%%%%%%%%%%%%%% STABILITIY CRITERIA %%%%%%%%%%%%%%%%%%%%%%%  

iter=iter+1
end



%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x=.5/nx:1/nx:1;
y=.5/ny:1/ny:1;
u(immersed) = nan;
v(immersed) = nan;
UX=reshape(u,nx,ny);
UY=reshape(v,nx,ny);
figure;streamslice(x,y,UX',UY')
title('Streams')
xlabel('x')
ylabel('y')
axis tight

[X,Y]=meshgrid(x,y);
[sx,sy] = meshgrid(.05,.05:.05:.95);
figure; 
streamline(stream2(X,Y,UX',UY',sx,sy))
title('Streamlines')
xlabel('x')
ylabel('y')

% Streamlines
figure;contourf(x,y,hypot(UX',UY'))

%%%%plot rho%%%%%%
rho(immersed) = nan;
rho_plot=transpose(squeeze(rho));
contourf(x, y, rho_plot, 'LineColor', 'none');
colorbar; % 显示颜色条
xlabel('x');
ylabel('y');
title('Rho Contour Graph');
figure;
%%%%plot pressure%%%%%%
pressure(immersed) = nan;
pressure_plot=transpose(squeeze(pressure));
contourf(x, y, pressure_plot, 'LineColor', 'none');
colorbar; % 显示颜色条
xlabel('x');
ylabel('y');
title('Pressure Contour Graph');
figure;

%%%%plot pressure%%%%%%
u(immersed) = nan;
u_plot=transpose(squeeze(u));
contourf(x, y, u_plot, 'LineColor', 'none');
colorbar; % 显示颜色条
xlabel('x');
ylabel('y');
title('u Contour Graph');
figure;
%%%%plot tau_alpha%%%%%%
tau_alpha(immersed) = nan;
tau_alpha_plot=transpose(squeeze(tau_alpha));
contourf(x, y, tau_alpha_plot, 'LineColor', 'none');
colorbar; % 显示颜色条
xlabel('x');
ylabel('y');
title('tau_alpha Contour Graph');

%%%%plot kne%%%%%%
kne(immersed) = nan;
kne_plot=transpose(squeeze(kne));
contourf(x, y, kne_plot, 'LineColor', 'none');
colorbar; % 显示颜色条
xlabel('x');
ylabel('y');
title('kne Contour Graph');
