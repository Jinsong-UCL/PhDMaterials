clear all

%load smile_fft.txt
%par=smile_fft;

load param_rec.txt
par=param_rec;


load param_recH0.txt
parH=param_recH0;


%nsize=max(size(fcallL(:,2)));


%% numero di osservazioni
nobs=36;

%% caso sei giorni
nobs=60;
%% numero strike
 nstr=5;
nv=12;
   
         maxdist=264-22
      
%% riordino dati per  scadenza e tipo
icont=0;
for io=1:nobs
 
      timetomat(io)=(maxdist-io+1);
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

figure(1)
h11=subplot(3,2,1)
plot(timetomat,epsr,'-','LineWidth',3)
ylabel('\eta ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat)) max(max(timetomat)) 0.15 0.251]) 


h22=subplot(3,2,2)
plot(timetomat,r0,'-','LineWidth',3)
ylabel('r_0  ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.007 0.0091])


h33=subplot(3,2,3)
plot(timetomat,chir,'-','LineWidth',3)
ylabel('\lambda   ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 3.0 4.5])

h44=subplot(3,2,4)
plot(timetomat,thetar,'-','LineWidth',3)
ylabel('\theta ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) 0 0.003])

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
axis([min(min(timetomat))  max(max(timetomat)) 0.8 1.2])


figure(2)
h1H=subplot(3,2,1)
plot(timetomat,epsv,'-',timetomat,epsvH,':','LineWidth',3)
ylabel('\gamma   ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat)) max(max(timetomat)) 0.7 1.0]) 
%pause(0.2)

h2H=subplot(3,2,2)
plot(timetomat,v0,'-',timetomat,v0H,':','LineWidth',3)
ylabel('v_0  ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.2 0.61])

h3H=subplot(3,2,3)
plot(timetomat,chiv,'-',timetomat,chivH,':','LineWidth',3)
ylabel('\chi   ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 5 9])

h4H=subplot(3,2,4)
plot(timetomat,thetav,'-',timetomat,thetavH,':','LineWidth',3)
ylabel('v^* ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) 0.00 0.02])

h5H=subplot(3,2,5)
plot(timetomat,rhov,'-',timetomat,rhovH,':','LineWidth',3)
ylabel('\rho_{p,v}     ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) -0.151 0.051])
%axis([min(min(timetomat))  max(max(timetomat)) -0.5 0.5])


h6H=subplot(3,2,6)
plot(timetomat,delta,'-',timetomat,deltaH,':','LineWidth',3)
ylabel('\Delta ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
%
axis([min(min(timetomat))  max(max(timetomat)) 0. 0.05])



