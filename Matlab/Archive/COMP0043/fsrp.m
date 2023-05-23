%% Monte Carlo simulation of the Feller square-root process
% dX = alpha*(mu-X)*dt + sigma*sqrt(X)*dW
% Used in the Cox-Ingersoll-Ross model and in the Heston stochastic volatility model

% Define parameters and time grid
clear variables % clear all variables from workspace
npaths = 2000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
alpha = 5; mu = 0.07; sigma = 0.265; % model parameters
%alpha = 5; mu = 0.03; sigma = 0.8; % model parameters
X0 = 0.03; % initial value
Feller_ratio = 2*alpha*mu/sigma^2 % for monitoring

%% Monte Carlo

% Allocate and initialise all paths
X = [X0*ones(1,npaths);zeros(nsteps,npaths)];

tic
% Euler-Maruyama
N = randn(nsteps,npaths); % sample standard normal random numbers
a = sigma^2/alpha*(exp(-alpha*dt)-exp(-2*alpha*dt)); % with analytic moments
b = mu*sigma^2/(2*alpha)*(1-exp(-alpha*dt))^2; % with analytic moments
for i = 1:nsteps % compute and accumulate the increments
   %X(i+1,:) = X(i,:) + alpha*(mu-X(i,:))*dt + sigma*sqrt(X(i,:)*dt).*N(i,:); % plain
    X(i+1,:) = mu+(X(i,:)-mu)*exp(-alpha*dt) + sqrt(a*X(i,:)+b).*N(i,:); % with analytic moments
   %if X(i+1,:) < 0
   %   X(i+1,:) = 0;
   %end
   %X(i+1,:) = X(i+1,:).*(X(i+1,:)>=0);
    X(i+1,:) = max(X(i+1,:),zeros(1,npaths));
end

% Exact method
% d = 4*alpha*mu/sigma^2; % degrees of freedom of the non-central chi square distribution
% k = sigma^2*(1-exp(-alpha*dt))/(4*alpha);
% for i = 1:nsteps % compute and accumulate the increments
%    lambda = 4*alpha*X(i,:)/(sigma^2*(exp(alpha*dt)-1));
%   %X(i+1,:) = icdf('ncx2',rand(1,npaths),d,lambda)*k; i % 80000 times slower than EM
%    X(i+1,:) = ncx2rnd(d,lambda,1,npaths)*k; % 40 times slower than EM
% end
toc

%% Expected, mean and sample paths
close all
figure(1)
EX = mu + (X0-mu)*exp(-alpha*t); % expected path
plot(t,EX,'k',t,mean(X,2),':k',t,mu*ones(size(t)),'k--',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path','\mu')
xlabel('t')
ylabel('X')
sdevinfty = sigma*sqrt(mu/(2*alpha));
ylim([-0.02 mu+4*sdevinfty])
title('Paths of a Feller square-root process dX = \alpha(\mu-X)dt + \sigmaX^{1/2}dW')
print('-dpng','fsrppaths.png')

%% Probability density function at different times
t2 = [0.05 0.1 0.2 0.4 1];
x = linspace(-0.02,mu+4*sdevinfty,200);
k = sigma^2*(1-exp(-alpha*t2))/(4*alpha);
d = 4*alpha*mu/sigma^2;
lambda = 4*alpha*X0./(sigma^2*(exp(alpha*t2)-1)); % non-centrality parameter
fa = zeros(length(x),length(t2)); % analytical
fs = zeros(length(x),length(t2)); % sampled
for i = 1:length(t2)
    fa(:,i) = pdf('ncx2',x/k(i),d,lambda(i))/k(i);
    fs(:,i) = hist(X(t2(i)*nsteps,:),x)/(npaths*(x(2)-x(1)));
end
figure(2)
plot(x,fa,x,fs)
xlabel('x')
ylabel('f_X(x,t)')
legend('t = 0.05','t = 0.10','t = 0.20','t = 0.40','t = 1.00')
title('Probability density function of a Feller square-root process at different times')
print('-dpng','fsrpdensities.png')
