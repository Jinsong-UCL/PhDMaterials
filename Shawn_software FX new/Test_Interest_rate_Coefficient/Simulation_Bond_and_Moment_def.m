clear all;
format long;
%%% four different interest rates
%% model m=0
fiduu=fopen('Result_moments_GRZ_1000_set_B_def10.txt','w'); 
fidu=fopen('Result_bond_GRZ_1000_set_B_def10.txt','w'); 
fidu1=fopen('Result_bond_GRZ_1000_set_B_def10_Table.txt','w'); 

randn('state',0);




%% Set A -  parameters
%% Parametri GRZ_OO_2011
% parameters of the volatility model
delta=0.01;
rho_v=-0.3;
chi_v=0.3;
theta_v=0.05;
eps_v=0.6;

%% parameters of the interest rate model
chi_r=0.01;
theta_r=0.02;
eps_r=0.01;
rho_r=-0.23
Omega=1.0;




%% Set A -  parameters
%% Parametri GRZ_OO_2011
% parameters of the volatility model
delta=0.01;
rho_v=-0.3;
chi_v=0.3;
theta_v=0.05;
eps_v=0.6;

%% parameters of the interest rate model
chi_r=0.01;
theta_r=0.02;
eps_r=0.01;
rho_r=-0.23
Omega=1.0;


   
%! param(1)=delta   		 
%! param(2) = eps_1 (vol of vol)
%! param (3)= v0_1 (volatility)
%! param(4)=theta_1
%!param(5)=chi_1
%! param (6)= rho_1 (correlation volatility - asset)
%! param(7)= eps_2 (vol of vol)
%! param(8) =v0_2 (volatility)
%! param(9)=theta_2
%!param(10)=chi_2
%! param(11)= rho_2 (correlation coefficient)

       
%% Set B -  parameters (coming from model calibration)
% parameters of the volatility model
delta=1.98;
rho_v=-0.97;
chi_v=0.65;
theta_v=0.0345;
eps_v=0.018;

%% parameters of the interest rate model
chi_r=19.37;
chi_r=3.62;
theta_r=0.00044
eps_r=0.0098;

rho_r=-0.81
Omega=2.51;

v0=0.089;
r0=0.00022;


2*theta_r*chi_r/eps_r^2
pause

n=4;

mtilde=1/2;
%%% Model like CIR


tmin=0;
tmax=0.25;
tmax=0.75;
tmax=2.0;

%%% Calcolo momenti uno e due per vari t

npoint=10;
%npoint=5;
tcurr(1)=0.25;
tcurr(2)=0.5;
tcurr(3)=0.75;
tcurr(4)=1.0
tcurr(5)=1.5;
tcurr(6)=2.0;
tcurr(7)=3.0;
tcurr(8)=5.0;
tcurr(9)=10.0;
tcurr(10)=20.0;


nsample=1000;
nsample=10;

 for i=1:npoint
    momuno(i)=0.0;
momdue(i)=0.0;

Bondtrue(i)=0.0;

end

for is=1:nsample
    
   is




for i=1: npoint

    
  
%% Attenzione se cambi i valori qui li devi cambiare anche sotto    
%%% Set A
   S0=100;
v0=0.05;
r0=0.02;
x0=0.0;


  
%% Set B (new calibro dati)
       S0=12.456;
x0=0.0;
v0=0.089;
r0=0.00022;


       %dt=dtl/chi_r;

       dt=0.01;
    npoint2=floor(tcurr(i)/dt);
    integ=0.0;

    for jj=1:npoint2
       rnd(1)=randn;
       rnd(2)=randn;
       rnd(3)=randn;
       rnd(4)=randn;

    
   x0=x0+(r0-0.5*(1+delta^2+2*delta*rho_v)*v0-0.5*Omega^2*r0)*dt+v0^0.5*(rnd(1)*(1.0-rho_v^2)^0.5+rnd(2)*rho_v)*dt^0.5+delta*v0^0.5*rnd(2)*dt^0.5...
       +Omega*r0^0.5*dt^0.5*(rnd(4)*(1.0-rho_r^2)^0.5+rnd(3)*rho_r);
    
   v0=v0+chi_v*(theta_v-v0)*dt+eps_v*v0^0.5*rnd(2)*dt^0.5;
    
  r0=r0+chi_r*(theta_r-r0)*dt+eps_r*r0^0.5*rnd(3)*dt^0.5;
   %r0=r0+(theta_r-r0)*dtl+(eps_r/chi_r^0.5)*r0^0.5*rnd(3)*dtl^0.5;
   
      
  integ=integ+r0*dt;  
   end
  
   
   Bondtrue(i)=Bondtrue(i)+1.0/exp(integ);
   momuno(i)=momuno(i)+(S0*exp(x0));
momdue(i)=momdue(i)+(S0*exp(x0))^2;

 
end

 pause(0.01) 
   %% closure sample
end

for i=1:npoint
     Bondtrue(i)=Bondtrue(i)/nsample
     momuno(i)=momuno(i)/nsample;
     momdue(i)=momdue(i)/nsample;
end    




%%%% calcolo teorico
%% Set B
    S0=12.456
v0=0.089;
r0=0.00022;
x0=0.0;



%%%%%%%%%%%% INSER DATA HERE

%%% Set A
%%% Set A
   S0=100;
v0=0.05;
r0=0.02;
x0=0.0;

%% Set B (new calibro dati set b def 1)
       S0=12.456;

x0=0.0;

v0=0.089;
r0=0.00022;

integ=0.;
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
  
 
nubd=2*chi_r*theta_r/eps_r^2;
srgbd=(1-1.0/exp(chi_r*tau));

%Mlambd=(4.0*chi_r/(4*chi_r+tau*eps_r^2*srgbd));
%Bondapp(i)=exp(-0.5*r0*tau)*Mlambd^nubd*exp(-0.5*Mlambd*tau*r0*exp(-chi_r*tau));


cc=eps_r^2*srgbd/(2.0*chi_r);

qq=(-exp(-chi_r*tau)+1);

qq=tau;
%Mlambd=(1.0/(1.0+tau*cc));

Mlambd=(1.0/(1.0+qq*cc));

%Bondapp(i)=1.0/exp(nubd*log(1.0/Mlambd)+Mlambd*r0*(exp(chi_r*tau)-1)/(chi_r*exp(chi_r*tau)));



%%% quarta formula
%nubd=2*chi_r*theta_r/eps_r^2;
%srgbd=(1-1.0/exp(chi_r*tau));
%au=tau^2*eps_r^2*0.25;
%cc=eps_r^2*srgbd/(2.0*chi_r);
%Mlambd=(1.0/(1.0+au*cc));
%Bondapp(i)=1.0/exp(nubd*log(1.0/Mlambd)+Mlambd*r0*au*exp(-chi_r*tau)+theta_r*tau+(r0-theta_r)*srgbd/chi_r);

%  uso media + integrazione

%% quinta formula
%nubd=2*chi_r*theta_r/eps_r^2;
%srgbd=(1-1.0/exp(chi_r*tau));
%au=(tau-srgbd/chi_r);
%cc=eps_r^2*srgbd/(2.0*chi_r);
%Mlambd=(1.0/(1.0+au*cc));
%Bondapp(i)=exp(-r0*srgbd/chi_r)*Mlambd^nubd*exp(-Mlambd*au*r0*(exp(-chi_r*tau)));


%% uso della media
%Bondappm(i)=1.0/exp(theta_r*tau+(r0-theta_r)*srgbd/chi_r);



%% sesta formula
nubd=2*chi_r*theta_r/eps_r^2;
srgbd=(1-1.0/exp(chi_r*tau));
au=tau*exp(chi_r*tau)/(1+exp(chi_r*tau));

hh=(2.0*eps_r^2+chi_r^2)^0.5;

au=tau*exp(0.5*hh*tau)/(1+exp(0.5*hh*tau));


cc=eps_r^2*srgbd/(2.0*chi_r);
Mlambd=(1.0/(1.0+au*cc));
Bondapp(i)=exp(-r0*tau/(1+exp(0.5*hh*tau)))*Mlambd^nubd*exp(-Mlambd*au*r0*(exp(-chi_r*tau)));

%Bondapp(i)=exp(-r0*tau/(1+chi_r))*Mlambd^nubd*exp(-Mlambd*au*r0*(exp(-chi_r*tau)));

%% uso della media
Bondappm(i)=1.0/exp(theta_r*tau+(r0-theta_r)*(srgbd/chi_r));



h=(chi_r^2+2.0*eps_r^2)^0.5;
aux=1.0-exp(-h*tau);

AA=2.0*h*exp((chi_r-h)*tau*0.5)/(2*h*exp(-h*tau)+(h+chi_r)*aux);
BB=2.0*aux/(2*h*exp(-h*tau)+(h+chi_r)*aux);
Bondtruev(i)=AA^(nubd)*exp(-r0*BB);

epp(i)=abs(Bondtruev(i)-Bondapp(i))/Bondtruev(i);
eppm(i)=abs(Bondtruev(i)-Bondappm(i))/Bondtruev(i);

end

  
  
for i=1:npoint
   fprintf(fiduu,'%e %e %e %e %e %e %e\n',tcurr(i),mom1exp(i),momuno(i),mom2exp(i),momdue(i),err1(i),err2(i));
   sq=sprintf('%e %e %e %e %e %e %e\n',tcurr(i),mom1exp(i),momuno(i),mom2exp(i),momdue(i),err1(i),err2(i));
  % disp(sq);
   
     fprintf(fidu,'%e %e %e %e %e %e %e\n',tcurr(i), Bondtrue(i),Bondtruev(i),Bondapp(i),epp(i),Bondappm(i),eppm(i));
     fprintf(fidu1,'%e %e %e %e \n',tcurr(i), Bondtruev(i),Bondapp(i),epp(i));
     
     sq1= fprintf('%e %e %e %e %e %e %e \n',tcurr(i), Bondtrue(i),Bondtruev(i),Bondapp(i),epp(i),Bondappm(i),eppm(i));
  disp(sq1);
end

fclose(fiduu);
fclose(fidu);
   fclose(fidu1);