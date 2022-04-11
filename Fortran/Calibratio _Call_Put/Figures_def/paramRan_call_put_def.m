clear all

%load smile_fft.txt
%par=smile_fft;

load param_rec_call_put_15_02_def3.txt
par=param_rec_call_put_15_02_def3;


load param_recH0_call_put.txt
parH=param_recH0_call_put;


load price_in_sample_call_put_15_02_def3.txt
pr=price_in_sample_call_put_15_02_def3;



%load price_in_sample_call_put_v1.txt
%pr=price_in_sample_call_put_v1;

%nsize=max(size(fcallL(:,2)));


%% caso sei giorni
nobs=60;
%% numero strike
 nstr=5;
nv=12;
   
         maxdist=264-22-60
      
%% riordino dati per  scadenza e tipo
icont=0;
for io=1:nobs
 
      timetomatt(io)=(maxdist-io+1)/264;
      timetomat(io)=io;
      %/252.0;
      %/263.d0;
    %  timetomat(io)=(maxdist-io+1)/252.d0;
      delta(io)=par(icont+1,2)
      epsv(io)=par(icont+2,2);
      v0(io)=par(icont+3,2);
      thetav(io)=par(icont+4,2);
      chiv(io)=par(icont+5,2);
      rhov(io)=par(icont+6,2);
      epsr(io)=par(icont+7,2);
      r0(io)=par(icont+8,2);
      thetar(io)=par(icont+9,2);
      chir(io)=par(icont+10,2);
      rhor(io)=par(icont+11,2);
     omega(io)=par(icont+12,2)
      
       deltaH(io)=parH(icont+1,2)
      epsvH(io)=parH(icont+2,2);
      v0H(io)=parH(icont+3,2);
      thetavH(io)=parH(icont+4,2);
      chivH(io)=parH(icont+5,2);
      rhovH(io)=parH(icont+6,2);
      epsrH(io)=parH(icont+7,2);
      r0H(io)=parH(icont+8,2);
      thetarH(io)=parH(icont+9,2);
      chirH(io)=parH(icont+10,2);
      rhorH(io)=parH(icont+11,2);
     omegaH(io)=parH(icont+12,2)
    
     
      icont=icont+nv;
end 

icont=0;
for io=1:nobs
    for j=1:6
        for i=1:5
        icont=icont+1
        prS0(i,io)=pr(icont,2);
        callT(i,io)=pr(icont,4);
         callA(i,io)=pr(icont,3);
          putT(i,io)=pr(icont,6);
         putA(i,io)=pr(icont,5);
     end
       %% nelle lettura resta memorizzato l'ultimo della finestra 
        
    end
end

figure(1)

h11=subplot(3,2,1)
plot(timetomat,epsr,'-','LineWidth',3)
ylabel('\eta ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat)) max(max(timetomat)) 0.05 0.1]) 


h22=subplot(3,2,2)
plot(timetomat,r0,'-',timetomat,r0H,':','LineWidth',3)
ylabel('r_0  ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
%axis([min(min(timetomat))  max(max(timetomat)) 0.8 1.21])


h33=subplot(3,2,3)
plot(timetomat,chir,'-','LineWidth',3)
ylabel('\lambda   ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 18 20])

h44=subplot(3,2,4)
plot(timetomat,thetar,'-','LineWidth',3)
ylabel('\theta ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.003 0.005])

h5=subplot(3,2,5)
plot(timetomat,rhor,'-','LineWidth',3)
ylabel('\rho_{p,r}     ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) -1.0 -0.8])

h6=subplot(3,2,6)
plot(timetomat,omega,'-','LineWidth',3)
ylabel('\Omega ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.5 1.5])



figure(2)

h1H=subplot(3,2,1)
plot(timetomat,epsv,'-',timetomat,epsvH,':','LineWidth',3)
ylabel('\gamma   ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat)) max(max(timetomat)) 0.0 0.31]) 
%pause(0.2)

h2H=subplot(3,2,2)
plot(timetomat,v0,'-',timetomat,v0H,':','LineWidth',3)
ylabel('v_0  ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.0 0.31])

h3H=subplot(3,2,3)
plot(timetomat,chiv,'-',timetomat,chivH,':','LineWidth',3)
ylabel('\chi   ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0 9])

h4H=subplot(3,2,4)
plot(timetomat,thetav,'-',timetomat,thetavH,':','LineWidth',3)
ylabel('v^* ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.00 0.31])

h5H=subplot(3,2,5)
plot(timetomat,rhov,'-',timetomat,rhovH,':','LineWidth',3)
ylabel('\rho_{p,v}     ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
%axis([min(min(timetomat))  max(max(timetomat)) -0.151 0.051])


h6H=subplot(3,2,6)
plot(timetomat,delta,'-',timetomat,deltaH,':','LineWidth',3)
ylabel('\Delta ','Rotation',0,'FontSize',16)
xlabel('    day index','FontSize',14)
%
axis([min(min(timetomat))  max(max(timetomat)) 0. 1.5])



figure(3)
ik=1
h1call=subplot(2,1,1)
plot(timetomat,callT(ik,:),'-',timetomat,callA(ik,:),':','LineWidth',3)
ylabel('C ','Rotation',0,'FontSize',16)
xlabel('    day index','FontSize',14)
legend('True call value', 'Approximated call value')
%
%axis([min(min(timetomat))  max(max(timetomat)) 0. 0.05])
h1put=subplot(2,1,2)
plot(timetomat,putT(ik,:),'-',timetomat,putA(ik,:),':','LineWidth',3)
ylabel('P ','Rotation',0,'FontSize',16)
xlabel('    day index','FontSize',14)
legend('True call value', 'Approximated call value')
%
%axis([min(min(timetomat))  max(max(timetomat)) 0. 0.05])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% confroto tassi
for i=1:max(size(timetomatt))
    mediar(i)=(mean(r0)-mean(thetar))*exp(-mean(chir)*(180-i)/264)+mean(thetar)*(1-exp(-mean(chir)*(180-i)/264));
    mediarH(i)=mean(r0H);
end
figure(4)
plot(timetomat,mediar,'-',timetomat,mediarH,':','LineWidth',3)
legend('Hybrid model expected value of interest rate','Heston model interest rate')

