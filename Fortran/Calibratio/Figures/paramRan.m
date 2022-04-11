clear all

%load smile_fft.txt
%par=smile_fft;

load smile_fft_ptime6days.txt
par=smile_fft_ptime6days;


%nsize=max(size(fcallL(:,2)));


%% numero di osservazioni
nobs=36;

%% caso sei giorni
nobs=30;
%% numero strike
 nstr=5;
nv=11;
   
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
     
      
      icont=icont+nv;
end 

figure(1)
h1=subplot(3,2,1)
plot(timetomat,epsv,'-s','LineWidth',2)
ylabel('\epsilon_v ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity (days)','FontSize',14)
xlabel('    day index','FontSize',14)

%xlabel('days from September 14th, 2010','FontSize',10)
%axis([208  238 0 0.5])

axis([min(min(timetomat)) max(max(timetomat)) 0.75 0.85]) 
%pause(0.2)

h2=subplot(3,2,2)
plot(timetomat,v0,'-s','LineWidth',2)
ylabel('v_0  ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
%xlabel('days from September 14th, 2010','FontSize',10)
axis([min(min(timetomat))  max(max(timetomat)) 0.35 0.65])

h3=subplot(3,2,3)
plot(timetomat,chiv,'-s','LineWidth',2)
ylabel('\chi_v ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) 3.5 6.5])

h4=subplot(3,2,4)
plot(timetomat,thetav,'-s','LineWidth',2)
ylabel('\theta_v ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) 0.008 0.012])

h5=subplot(3,2,5)
plot(timetomat,rhov,'-s','LineWidth',2)
ylabel('\rho_v ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) -0.5 0.5])


h6=subplot(3,2,6)
plot(timetomat,delta,'-s','LineWidth',2)
ylabel('\Delta ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) 0.4 1.0])




figure(2)
h11=subplot(3,2,1)
plot(timetomat,epsr,'-s','LineWidth',2)
ylabel('\epsilon_r ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat)) max(max(timetomat)) 0.055 0.065]) 


h22=subplot(3,2,2)
plot(timetomat,r0,'-s','LineWidth',2)
ylabel('r_0  ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 0.0055 0.0095])


h33=subplot(3,2,3)
plot(timetomat,chir,'-s','LineWidth',2)
ylabel('\chi_r ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)
axis([min(min(timetomat))  max(max(timetomat)) 3.5 6.5])

h44=subplot(3,2,4)
plot(timetomat,thetar,'-s','LineWidth',2)
ylabel('\theta_r ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) 0 0.04])

h5=subplot(3,2,5)
plot(timetomat,rhor,'-s','LineWidth',2)
ylabel('\rho_r ','Rotation',0,'FontSize',16)
%xlabel('    time to maturity  (days)','FontSize',14)
xlabel('    day index','FontSize',14)

axis([min(min(timetomat))  max(max(timetomat)) -1 -0.8])







