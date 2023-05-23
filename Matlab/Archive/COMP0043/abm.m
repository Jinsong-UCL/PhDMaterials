%% Monte Carlo simulation of arithmetic Brownian motion
% dX(t) = mu*dt + sigma*dW(t)

% Define parameters and time grid
clear variables % clear all variables from workspace
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
%mu = 0.12; sigma = 0.4; % model parameters
mu = -0.05; sigma = 0.4; % model parameters

%% Monte Carlo

% Compute the increments with Euler-Maruyama
dX = mu*dt + sigma*randn(npaths,nsteps)*sqrt(dt);

% Accumulate the increments
X = [zeros(npaths,1) cumsum(dX,2)];

%% Expected, mean and sample path
close all
figure(1)
EX = mu*t; % expected path
plot(t,EX,'k',t,mean(X),':k',t,X(1:1000:end,:),t,EX,'k',t,mean(X),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
ylim([-1,1]);
title('Arithmetic Brownian motion dX(t) = \mudt + \sigmadW(t)')
print('-dpng','abmpaths.png')

%% Variance = mean square deviation = mean square displacement of the random part
figure(2)
plot(t,sigma^2*t,t,var(X))
legend('Theory: \sigma^2t = 2Dt','Sampled','Location','NorthWest')
xlabel('t')
ylabel('Var(X) = E((X-E(X))^2)')
title('Arithmetic Brownian motion: MSD')
print('-dpng','abmmsd.png')

%% Mean absolute deviation
figure(3)
plot(t,sigma*sqrt(2*t/pi),t,mean(abs(X-EX)))
legend('Theory: \sigma(2t/\pi)^{1/2}','Sampled','Location','NorthWest')
xlabel('t')
ylabel('E(|X-E(X)|) = (2Var(X)/pi)^{1/2}')
ylim([0 0.02])
title('Arithmetic Brownian motion: mean absolute deviation')
print('-dpng','mad.png')

%% Probability density function at different times
figure(4)

subplot(3,1,1)
histogram(X(:,20),-1:0.02:1,'normalization','pdf');
ylabel('f_X(x,0.1)')
xlim([-1,1])
ylim([0,3.5])
title('Arithmetic Brownian motion: PDF at different times')

subplot(3,1,2)
histogram(X(:,80),-1:0.02:1,'normalization','pdf');
xlim([-1,1])
ylim([0,3.5])
ylabel('f_X(x,0.4)')

subplot(3,1,3)
histogram(X(:,end),-1:0.02:1,'normalization','pdf');
xlim([-1,1])
ylim([0,3.5])
xlabel('x')
ylabel('f_X(x,1)')
print('-dpng','abmhist.png')

%% Solution of the Fokker-Planck equation
figure(5)
D = sigma^2/2; % diffusion coefficient
[x,tt] = meshgrid(-1:0.02:1,0.1:0.025:1);
f = 1./(2*sqrt(pi*D*tt)).*exp(-(x-mu*tt).^2./(4*D*tt));
mesh(x,tt,f)
%contour(x,tt,f,100)
xlabel('x')
ylabel('t')
zlabel('f_X')
title('Arithmetic Brownian motion: solution of the Fokker-Planck equation with \mu = -0.05, \sigma = 0.4')
view(30,24)
print('-dpng','abmfpe.png')
