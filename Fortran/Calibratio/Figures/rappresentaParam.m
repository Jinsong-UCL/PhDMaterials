clear all

load Forecast_call_00.txt
load Forecast_put_00.txt

fcallL=Forecast_call_00;
fputL=Forecast_put_00;

futS0=fcallL(:,1);
tmat=fcallL(:,2);
strike=fcallL(:,3);

forecall=fcallL(:,4);
obscall=fcallL(:,5);

foreput=fputL(:,4);
obsput=fputL(:,5);


nsize=max(size(fcallL(:,2)));


%% numero di osservazioni
nobs=30;
%% numero strike
 nstr=5;
%
%s0=1.3446;
%s1=1.3536

%% riordino dati per strike e scadenza
icont=0;
for io=1:nobs
  for is=1:nstr
    icont=icont+1;
      timetomat(is,io)=tmat(icont)*252;
      mstrike(is,io)=strike(icont)/futS0(icont);
      callobs(is,io)=obscall(icont);
      callfor(is,io)=forecall(icont);
      errcall(is,io)=abs(obscall(icont)-forecall(icont))/abs(obscall(icont));
%% ora procedo per la put
   putobs(is,io)=obsput(icont);
      putfor(is,io)=foreput(icont);
      errput(is,io)=abs(obsput(icont)-foreput(icont))/abs(obsput(icont));

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
axis([208 235 10 400])

subplot(3,2,2)
is=1;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(a)    time to maturity (days)','FontSize',16)
legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([208 235 10 400])

subplot(3,2,3)
is=2;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(b)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([208 235 10 400])

subplot(3,2,4)
is=3;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(c)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([208 235 10 400])

subplot(3,2,5)
is=4;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(d)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([208 235 10 400])

subplot(3,2,6)
is=5;
plot(timetomat(is,n1:n2),callobs(is,n1:n2),'-s',timetomat(is,n1:n2),callfor(is,n1:n2),':p',timetomat(is,n1:n2),putobs(is,n1:n2),'-s',timetomat(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(e)    time to maturity (days)','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
axis([208 235 10 400])


figure(2)

subplot(2,2,1)
is=1;
plot(mstrike(is,n1:n2),callobs(is,n1:n2),'-s',mstrike(is,n1:n2),callfor(is,n1:n2),':p',mstrike(is,n1:n2),putobs(is,n1:n2),'-s',mstrike(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(a)  moneyness ','FontSize',16)
legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
%axis([208 235 10 350])

subplot(2,2,2)
is=2;
plot(mstrike(is,n1:n2),callobs(is,n1:n2),'-s',mstrike(is,n1:n2),callfor(is,n1:n2),':p',mstrike(is,n1:n2),putobs(is,n1:n2),'-s',mstrike(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(b)    moneyness','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
%axis([208 235 10 350])

subplot(2,2,3)
is=3;
plot(mstrike(is,n1:n2),callobs(is,n1:n2),'-s',mstrike(is,n1:n2),callfor(is,n1:n2),':p',mstrike(is,n1:n2),putobs(is,n1:n2),'-s',mstrike(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(c)    moneyness','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
%axis([208 235 10 350])

subplot(2,2,4)
is=4;
plot(mstrike(is,n1:n2),callobs(is,n1:n2),'-s',mstrike(is,n1:n2),callfor(is,n1:n2),':p',mstrike(is,n1:n2),putobs(is,n1:n2),'-s',mstrike(is,n1:n2),putfor(is,n1:n2),':p','LineWidth',2)
ylabel('V  ','Rotation',0,'FontSize',16)
xlabel('(d)   moneyness','FontSize',16)
%legend('Observed call option prices', 'Forecast call option prices','Observed put option prices', 'Forecast put option prices')
%axis([208 235 10 350])