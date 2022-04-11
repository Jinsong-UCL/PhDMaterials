clear all;
%%% four different interest rates
%% model m=0
%fiduu=fopen('Result_moments.txt','w'); 

%load Result_moments_ottimo.txt
%rr=Result_moments_ottimo;

load Result_moments_GRZ_1000.txt
rr=Result_moments_GRZ_1000;


tcurr=rr(:,1);
mom_one_app=rr(:,3);
mom_one_true=rr(:,2);

mom_two_app=rr(:,5);
mom_two_true=rr(:,4);
err_m1=rr(:,6);
err_m2=rr(:,7);


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


%% Attenzione se cambi i valori qui li devi cambiare anche sotto    
    S0=12.456;
v0=0.25;
r0=0.13;
x0=0.0;


figure(1)
h1=subplot(2,2,1)
plot(tcurr,mom_one_true,'-',tcurr,mom_one_app,':')
xlabel('time')
ylabel('M_1    ','Rotation',0)
legend('True M_1','Monte Carlo M_1',2)
%axis([0 1 12 16])
axis([0 1 96 116])


h2=subplot(2,2,2)
plot(tcurr,mom_two_true,'-',tcurr,mom_two_app,':')
xlabel('time')
ylabel('M_2    ','Rotation',0)
legend('True M_2','Monte Carlo M_2',2)
axis([0 1 9500 13000])
%axis([0 1 160 280])
%legend('True second order moment','MC second order moment')


h3=subplot(2,2,3)
plot(tcurr,err_m1,'-')
xlabel('time')
ylabel('Relative Error M_1','Rotation',90)
axis([0 1 0 0.0301])
%legend('Relative error first order moment')

h4=subplot(2,2,4)
plot(tcurr,err_m2,'-')
xlabel('time')
ylabel('Relative Error M_2','Rotation',90)
axis([0 1 0 0.0301])
%legend('Relative error second order moment')

%for i=1:npoint
%   fprintf(fiduu,'%e %e %e %e %e %e %e\n',tcurr(i),mom1exp(i),momuno(i),mom2exp(i),momdue(i),err1(i),err2(i));
%   sq=sprintf('%e %e %e %e %e %e %e\n',tcurr(i),mom1exp(i),momuno(i),mom2exp(i),momdue(i),err1(i),err2(i));
%   disp(sq);
%end
%
%fclose(fiduu);
   