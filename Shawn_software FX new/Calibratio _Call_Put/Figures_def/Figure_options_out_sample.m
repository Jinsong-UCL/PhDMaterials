clear all


load out_sample_call_put_67_84_def3.txt
dati_all=out_sample_call_put_67_84_def3;


load out_sampleH0_call_put_67_84.txt
dati_allH0=out_sampleH0_call_put_67_84;

nptime=1;
nobs=17;



strike(1)=1075;
strike(2)=1100;
strike(3)=1125;
strike(4)=1150;
strike(5)=1170;


futS0=dati_all(:,2);

forecall=dati_all(:,3);
obscall=dati_all(:,4);

foreput=dati_all(:,5);
obsput=dati_all(:,6);


forecallH0=dati_allH0(:,3);
foreputH0=dati_allH0(:,5);



nsize=max(size(futS0));


%% numero di osservazioni
%% numero strike
 nstr=5;
nv=6;

  
	   maxdist=264-22-66;
%% riordino dati per strike e scadenza
icont=0;

for io=1:nobs
    
    tt(io)=(maxdist-nptime-io+1);
  
    if(nptime>1) 
      for it=1:nptime-1
      for is=1:nstr
        icont=icont+1;
      end
      end
  end
  
for is=1:nstr 
        icont=icont+1;

           mstrike(is,io)=strike(is)/futS0(icont);
           timetomat(is,io)=(maxdist-io+1);
      
      callobs(is,io)=obscall(icont);
      callfor(is,io)=forecall(icont);
      errcall(is,io)=abs(obscall(icont)-forecall(icont))/abs(obscall(icont));
      
   
      callforH0(is,io)=forecallH0(icont);
      errcallH0(is,io)=abs(obscall(icont)-forecallH0(icont))/abs(obscall(icont));
%% ora procedo per la put
   putobs(is,io)=obsput(icont);
      putfor(is,io)=foreput(icont);
      errput(is,io)=abs(obsput(icont)-foreput(icont))/abs(obsput(icont));
      
          putforH0(is,io)=foreputH0(icont);
      errputH0(is,io)=abs(obsput(icont)-foreputH0(icont))/abs(obsput(icont));
  %         icont=icont+1;      
   end
   
end 

%% mean error

    
%% 21 days is a month
n1=1;
n2=nobs;
for is=1:nstr
  ermeancall(is)=mean(errcall(is,n1:n2));
  ermeanput(is)=mean(errput(is,n1:n2));
end



disp('errore medio call one factor') 
mn(1)=mean(mean(errcall))
disp('errore medio put one factor')
mn(2)=mean(mean(errput))

pause

disp('errore medio call max factors') 
sum(sum(errcall))
disp('errore medio put max factors')
sum(sum(errput))

disp('errore medio call max factors') 
max(max(errcall))
disp('errore medio put max factors')
max(max(errput))



figure(1)

subplot(3,2,3)
is=2;
%plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
%ylabel('V  ','Rotation',0,'FontSize',16)
%xlabel('(a)    time to maturity (days)','FontSize',16)
legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices','FontSize',14)
axis([183 242 10 400])

subplot(3,2,2)
is=1;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(a)    time to maturity (days)','FontSize',16)
legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([183 242 10 400])

subplot(3,2,3)
is=2;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(b)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([183 242 10 400])

subplot(3,2,4)
is=3;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(c)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([183 242 10 400])

subplot(3,2,5)
is=4;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(d)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([183 242 10 400])

subplot(3,2,6)
is=5;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(e)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([183 242 10 400])


init=160;
inif=176;

figure(2)


subplot(3,2,3)
is=2;
%plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
%ylabel('V  ','Rotation',0,'FontSize',16)
%xlabel('(a)    time to maturity (days)','FontSize',16)
legend('Observed call option prices', 'Hybrid model forecast call option prices', 'Heston model forecast call option prices','FontSize',14)
axis([init inif 150 455])

subplot(3,2,2)
is=1;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-',timetomat(is,n1:n2),callfor(is,n1:n2),':',timetomat(is,n1:n2),callforH0(is,n1:n2),'--','LineWidth',2)
ylabel('C  ','Rotation',0,'FontSize',16)
xlabel('(a)    time to maturity (days)','FontSize',16)
legend('Observed call option prices', 'Hybrid model forecast call option prices', 'Heston model forecast call option prices')
axis([init inif 150 455])

subplot(3,2,3)
is=2;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-',timetomat(is,n1:n2),callfor(is,n1:n2),':',timetomat(is,n1:n2),callforH0(is,n1:n2),'--','LineWidth',2)
ylabel('C  ','Rotation',0,'FontSize',16)
xlabel('(b)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([init inif 150 455])

subplot(3,2,4)
is=3;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-',timetomat(is,n1:n2),callfor(is,n1:n2),':',timetomat(is,n1:n2),callforH0(is,n1:n2),'--','LineWidth',2)
%
%plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('C  ','Rotation',0,'FontSize',16)
xlabel('(c)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([init inif 150 455])

subplot(3,2,5)
is=4;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-',timetomat(is,n1:n2),callfor(is,n1:n2),':',timetomat(is,n1:n2),callforH0(is,n1:n2),'--','LineWidth',2)
%
%plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('C  ','Rotation',0,'FontSize',16)
xlabel('(d)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([init inif 150 455])

subplot(3,2,6)
is=5;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-',timetomat(is,n1:n2),callfor(is,n1:n2),':',timetomat(is,n1:n2),callforH0(is,n1:n2),'--','LineWidth',2)
%plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('C  ','Rotation',0,'FontSize',16)
xlabel('(e)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([init inif 150 455])





figure(3)


subplot(3,2,3)
is=2;
legend('Observed put option prices', 'Hybrid model forecast put option prices', 'Heston model forecast put option prices','FontSize',14)
axis([init inif 10 80])

subplot(3,2,2)
is=1;
plot(timetomat(is,n1:n2),putobs(is,n1:n2),'-',timetomat(is,n1:n2),putfor(is,n1:n2),':',timetomat(is,n1:n2),putforH0(is,n1:n2),'--','LineWidth',2)
ylabel('P  ','Rotation',0,'FontSize',16)
xlabel('(a)    time to maturity (days)','FontSize',16)
legend('Observed put option prices', 'Hybrid model forecast put option prices', 'Heston model forecast put option prices')
axis([init inif 10 80])

subplot(3,2,3)
is=2;
plot(timetomat(is,n1:n2),putobs(is,n1:n2),'-',timetomat(is,n1:n2),putfor(is,n1:n2),':',timetomat(is,n1:n2),putforH0(is,n1:n2),'--','LineWidth',2)
ylabel('P  ','Rotation',0,'FontSize',16)
xlabel('(b)    time to maturity (days)','FontSize',16)
axis([init inif 10 80])

subplot(3,2,4)
is=3;
plot(timetomat(is,n1:n2),putobs(is,n1:n2),'-',timetomat(is,n1:n2),putfor(is,n1:n2),':',timetomat(is,n1:n2),putforH0(is,n1:n2),'--','LineWidth',2)
ylabel('P  ','Rotation',0,'FontSize',16)
xlabel('(c)    time to maturity (days)','FontSize',16)
axis([init inif 10 80])

subplot(3,2,5)
is=4;
plot(timetomat(is,n1:n2),putobs(is,n1:n2),'-',timetomat(is,n1:n2),putfor(is,n1:n2),':',timetomat(is,n1:n2),putforH0(is,n1:n2),'--','LineWidth',2)
ylabel('P  ','Rotation',0,'FontSize',16)
xlabel('(d)    time to maturity (days)','FontSize',16)
axis([init inif 10 80])

subplot(3,2,6)
is=5;
plot(timetomat(is,n1:n2),putobs(is,n1:n2),'-',timetomat(is,n1:n2),putfor(is,n1:n2),':',timetomat(is,n1:n2),putforH0(is,n1:n2),'--','LineWidth',2)
ylabel('P  ','Rotation',0,'FontSize',16)
xlabel('(e)    time to maturity (days)','FontSize',16)
axis([init inif 10 80])

