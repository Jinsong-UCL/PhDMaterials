clear all;
%%% four different interest rates
%% model m=0
fiduu=fopen('Result_moments.txt','w'); 

randn('state',0);

% parameters of the volatility model
delta=0.01;
rho_v=-0.02;
chi_v=2.1;
theta_v=0.01;
eps_v=0.01;

%% parameters of the interest rate model
chi_r=2.1;
theta_r=0.01;
eps_r=0.01;
rho_r=-0.23
Omega=1.0;

% chi_r>(2-Omega)eta_r


n=4;

mtilde=1/2;
%%% Model like CIR


tmin=0;
tmax=0.25;
tmax=0.75;
tmax=1.0;

%%% Calcolo momenti uno e due per vari t

npoint=20000;
npoint=10;
npoint2=2000;
% time step
dt=(tmax-tmin)/(npoint*npoint2);


nsample=100;


Bondtrue=0.0;

 for i=1:npoint
    momuno(i)=0.0;
momdue(i)=0.0;
end

for is=1:nsample
    
   is

%% Attenzione se cambi i valori qui li devi cambiare anche sotto    
    S0=12.456;
v0=0.25;
r0=0.13;
x0=0.0;
integ=0.0;



dt=0.00001;

for i=1: npoint
    
       S0=12.456;
v0=0.25;
r0=0.13;
x0=0.0;
integ=0.0;
    
    tcurr(i)=i*tmax/npoint;

    npoint2=floor(tcurr(i)/dt);
    
    for jj=1:npoint2
       rnd(1)=randn;
       rnd(2)=randn;
       rnd(3)=randn;
       rnd(4)=randn;

    
   x0=x0+(r0-0.5*(1+delta^2+2*delta*rho_v)*v0-0.5*Omega^2*r0)*dt+v0^0.5*(rnd(1)*(1.0-rho_v^2)^0.5+rnd(2)*rho_v)*dt^0.5+delta*v0^0.5*rnd(2)*dt^0.5...
       +Omega*r0^0.5*dt^0.5*(rnd(4)*(1.0-rho_r^2)^0.5+rnd(3)*rho_r);
    
   v0=v0+chi_v*(theta_v-v0)*dt+eps_v*v0^0.5*rnd(2)*dt^0.5;
    
   r0=r0+chi_r*(theta_r-r0)*dt+eps_r*r0^0.5*rnd(3)*dt^0.5;
   
  
   end
  

   momuno(i)=momuno(i)+(S0*exp(x0));
momdue(i)=momdue(i)+(S0*exp(x0))^2;

 
end

 pause(0.01) 
   %% closure sample
end

for i=1:npoint
     momuno(i)=momuno(i)/nsample;
     momdue(i)=momdue(i)/nsample;
end    




%%%% calcolo teorico

    S0=12.456
v0=0.25;
r0=0.13;
x0=0.0;


mm=1.0;
muv1=-0.5*(chi_v-mm*eps_v*(delta+rho_v));
zetav1=0.5*(4.0*muv1^2-eps_v^2*(mm^2-mm)*(1+delta^2+2.0*delta*rho_v))^0.5;

mur1=-0.5*(chi_r-mm*eps_r*Omega*rho_r);
zetar1=0.5*(4.0*mur1^2-eps_r^2*((mm^2-mm)*Omega+2*mm))^0.5;


mm=2.0;
muv2=-0.5*(chi_v-mm*eps_v*(delta+rho_v));
zetav2=0.5*(4.0*muv2^2-eps_v^2*(mm^2-mm)*(1+delta^2+2.0*delta*rho_v))^0.5;

mur2=-0.5*(chi_r-mm*eps_r*Omega*rho_r);
zetar2=0.5*(4.0*mur2^2-eps_r^2*((mm^2-mm)*Omega+2*mm))^0.5;


for i=1:npoint
    tau=tcurr(i);
    
    
      zetav=zetav1;
      muv=muv1;
      zetar=zetar1;
      mur=mur1;
     
     sgv=1.0-exp(-2.0*tau*zetav);
     sgr=1.0-exp(-2.0*tau*zetar);

      sbv=(zetav+muv)*exp(-2.0*tau*zetav)+(zetav-muv);
      sbr=(zetar+mur)*exp(-2.0*tau*zetar)+(zetar-mur);
      
  p1v=2.0*chi_v*theta_v*log(sbv/(2.0*zetav))/eps_v^2;
  p2v=2.0*chi_v*theta_v*tau*(muv+zetav)/eps_v^2;
  p3v=2.0*v0*(zetav^2-muv^2)*sgv/(eps_v^2*sbv);
 
  
   p1r=2.0*chi_r*theta_r*log(sbr/(2.0*zetar))/eps_r^2;
  p2r=2.0*chi_r*theta_r*tau*(mur+zetar)/eps_r^2;
  p3r=2.0*r0*(zetar^2-mur^2)*sgr/(eps_r^2*sbr);
 
  mom1exp(i)=S0/exp(p1v+p2v+p3v+p1r+p2r+p3r);
  
    %%%%%%%%%%%5 compute second moment
    zetav=zetav2;
      muv=muv2;
      zetar=zetar2;
      mur=mur2;
     
     sgv=1.0-exp(-2.0*tau*zetav);
     sgr=1.0-exp(-2.0*tau*zetar);

      sbv=(zetav+muv)*exp(-2.0*tau*zetav)+(zetav-muv);
      sbr=(zetar+mur)*exp(-2.0*tau*zetar)+(zetar-mur);
      
  p1v=2.0*chi_v*theta_v*log(sbv/(2.0*zetav))/eps_v^2;
  p2v=2.0*chi_v*theta_v*tau*(muv+zetav)/eps_v^2;
  p3v=2.0*v0*(zetav^2-muv^2)*sgv/(eps_v^2*sbv);
 
  
   p1r=2.0*chi_r*theta_r*log(sbr/(2.0*zetar))/eps_r^2;
  p2r=2.0*chi_r*theta_r*tau*(mur+zetar)/eps_r^2;
  p3r=2.0*r0*(zetar^2-mur^2)*sgr/(eps_r^2*sbr);
 
  mom2exp(i)=S0^2/exp(p1v+p2v+p3v+p1r+p2r+p3r);

  
  err1(i)=abs(mom1exp(i)-momuno(i))/abs(momuno(i));
  err2(i)=abs(mom2exp(i)-momdue(i))/abs(momdue(i));
  
  end
  

for i=1:npoint
   fprintf(fiduu,'%e %e %e %e %e %e %e\n',tcurr(i),mom1exp(i),momuno(i),mom2exp(i),momdue(i),err1(i),err2(i));
   sq=sprintf('%e %e %e %e %e %e %e\n',tcurr(i),mom1exp(i),momuno(i),mom2exp(i),momdue(i),err1(i),err2(i));
   disp(sq);
end

fclose(fiduu);
   