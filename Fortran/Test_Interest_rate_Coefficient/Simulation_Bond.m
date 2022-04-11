clear all;
%%% four different interest rates
%% model m=0
fiduu=fopen('Result_app.txt','w'); 

randn('state',0);

% parameters of the volatility model
delta=0.01;
rho_v=0.02;
chi_v=2.1;
theta_v=0.01;
eps_v=0.01;

%% parameters of the interest rate model
chi_r=2.1;
theta_r=0.01;
eps_r=0.01;
rho_r=-0.23


n=4;

mtilde=1/2;
%%% Model like CIR


tmin=0;

tmax=0.75;
tmax=0.50;
%tmax=0.25;
%tmax=0.125;
%tmax=1.0/12.0;
tmax=1.0
npoint=40000;
%tmax=0.5
%npoint=10000;
%tmax=1.0/12.0;

%tmax=0.750
%npoint=20000;

% time step
dt=(tmax-tmin)/npoint;


nsample=1000;


Bondtrue=0.0;
for i=1:nsample

    
    S0=12.456
v0=0.25;
r0=0.13;
x0=0.0;
integ=0.0
for i=1: npoint
    
    rnd(1)=randn;
    rnd(2)=randn;
    rnd(3)=randn;
    rnd(4)=randn;
    
    
   x0=x0+(r0-0.5*(1+delta^2+2*delta*rho_v)*v0+r0)*dt+v0^0.5*(rnd(1)*(1.0-rho_v^2)^0.5+rnd(2)*rho_v)*dt^0.5+delta*rnd(2)*dt^0.5+r0^0.5*dt^0.5*(rnd(3)*(1.0-rho_r^2)^0.5+rnd(4)*rho_r);
    
   v0=v0+chi_v*(theta_v-v0)*dt+eps_v*v0^0.5*rnd(2)*dt^0.5;
    
   r0=r0+chi_r*(theta_r-r0)*dt+eps_r*r0^0.5*rnd(3)*dt^0.5;
   
   integ=integ+r0*dt;
   end
  
   
   Bondtrue=Bondtrue+1.0/exp(integ)
  pause(0.1) 
   %% closure sample
end
Bondtrue=Bondtrue/nsample;


    
    S0=12.456
v0=0.25;
r0=0.13;
x0=0.0;
nu=2*chi_r*theta_r/eps_r^2;
srg=(1-exp(-chi_r*tmax));
Mlam=(4.0*chi_r/(4*chi_r+tmax*eps_r^2*srg));
Bondapp=exp(-0.5*r0*tmax)*Mlam^nu*exp(-0.5*Mlam*tmax*r0*exp(-chi_r*tmax));

epp=abs(Bondtrue-Bondapp)/Bondtrue


   fprintf(fiduu,'%e %e %e %e %e %e \n',tmax,r0,v0, Bondtrue,Bondapp,epp);
fclose(fiduu);


   