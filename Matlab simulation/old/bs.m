%% Compute the Black-Scholes-Merton price of European options
% The model for the underlying is geometric Brownian motion
% dS = mu*S*dt + sigma*S*dW

% Contract parameters
T = 1; % maturity
K = 1.1; % strike price

% Market parameters
S0 = 1; % spot price
r = 0.05; % risk-free interest rate
q = 0.02; % dividend rate

% Model parameter
sigma = 0.4; % volatility

% Monte Carlo parameters; npaths = nblocks*nsample
nblocks = 2000; % number of blocks
nsample = 10000; % number of samples per block

% Fourier parameters
xwidth = 6; % width of the support in real space
ngrid = 2^8; % number of grid points
alpha = -10; % damping factor for a call

% Controls
figures = 0;

%% Analytical solution
tic
muABM = r-q-0.5*sigma^2; % drift coefficient of the arithmetic Brownian motion
d2 = (log(S0/K)+muABM*T)/(sigma*sqrt(T));
d1 = d2 + sigma*sqrt(T);
Vca = S0*exp(-q*T)*cdf('Normal',d1,0,1) - K*exp(-r*T)*cdf('Normal',d2,0,1);
Vpa = K*exp(-r*T)*cdf('Normal',-d2,0,1) - S0*exp(-q*T)*cdf('Normal',-d1,0,1);
% Put-call parity: Vp = Vc + K*exp(-(r-q)*T) + S0
cputime_a = toc;

% Analytical solution provided by Matlab's Financial Toolbox
tic
[VcaM,VpaM] = blsprice(S0,K,r,T,sigma,q);
cputime_aM = toc;

% Print the results
fprintf('%20s%14s%14s%14s\n','','call','put','CPU_time/s')
fprintf('%20s%14.10f%14.10f%14.10f\n','BS analytical',Vca,Vpa,cputime_a)
fprintf('%20s%14.10f%14.10f%14.10f\n','BS analytical Matlab',VcaM,VpaM,cputime_aM)

if figures ~= 0

    % Plot the analytical solution
    [St,t] = meshgrid(0:.05:2,0:0.025:T);
    d2 = (log(St/K)+muABM*(T-t))./(sigma*sqrt(T-t));
    d1 = d2 + sigma*sqrt(T-t);

    close all
    figure(1)
    Vc = St.*exp(-q*(T-t)).*cdf('Normal',d1,0,1) - K*exp(-r*(T-t)).*cdf('Normal',d2,0,1);
    Vc(end,:) = max(St(end,:)-K,0);
    mesh(St,t,Vc)
    xlabel('S')
    ylabel('t')
    zlabel('V')
    title('Call')
    view(-30,24)
    print('-dpng','bsc.png')
    
    figure(2)
    Vp = K*exp(-r*(T-t)).*cdf('Normal',-d2,0,1) - St.*exp(-q*(T-t)).*cdf('Normal',-d1,0,1);
    Vp(end,:) = max(K-St(end,:),0);
    mesh(St,t,Vp)
    xlabel('S')
    ylabel('t')
    zlabel('V')
    title('Put')
    view(30,24)
    print('-dpng','bsp.png')
    
    % Plot the analytical solution as a function of the log price
    k = log(K/S0);
    [xt,t] = meshgrid(-1:.05:1,0:0.025:T);
    d2 = (xt-k+muABM*(T-t))./(sigma*sqrt(T-t));
    d1 = d2 + sigma*sqrt(T-t);
    
    figure(3)
    Vc = S0*(exp(xt-q*(T-t)).*cdf('Normal',d1,0,1) - exp(k-r*(T-t)).*cdf('Normal',d2,0,1));
    Vc(end,:) = S0*max(exp(xt(end,:))-exp(k),0);
    mesh(xt,t,Vc)
    xlabel('x')
    ylabel('t')
    zlabel('V')
    title('Call')
    view(-30,24)
    print('-dpng','bscx.png')
    
    figure(4)
    Vp = S0*(exp(k-r*(T-t)).*cdf('Normal',-d2,0,1) - exp(xt-q*(T-t)).*cdf('Normal',-d1,0,1));
    Vp(end,:) = S0*max(exp(k)-exp(xt(end,:)),0);
    mesh(xt,t,Vp)
    xlabel('x')
    ylabel('t')
    zlabel('V')
    title('Put')
    view(30,24)
    print('-dpng','bspx.png')

end

%% Fourier transform method

% Grids in real and Fourier space
tic
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

% Characteristic function at time T
xia = xi+1i*alpha; % call
psi = 1i*muABM*xia-0.5*(sigma*xia).^2; % characteristic exponent
Psic = exp(psi*T); % characteristic function
xia = xi-1i*alpha; % put
psi = 1i*muABM*xia-0.5*(sigma*xia).^2; % characteristic exponent
Psip = exp(psi*T); % characteristic function

% These functions provide the characteristic functions of 8 Levy processes
% param = parameters(1,T,T,r,q); % set the parameters editing parameters.m
% [x,fc,xi,Psic] = kernel(ngrid,-b,b,param,alpha,0,1); % call
% [x,fp,xi,Psip] = kernel(ngrid,-b,b,param,-alpha,0,1); % put

% Fourier transform of the payoff
b = xwidth/2; % upper bound of the support in real space
U = S0*exp(b);
L = S0*exp(-b);
[~,gc,Gc] = payoff(x,xi,alpha,K,L,U,S0,1); % call
[S,gp,Gp] = payoff(x,xi,-alpha,K,L,U,S0,0); % put

% Discounted expected payoff computed with the Plancherel theorem
c = exp(-r*T).*real(fftshift(fft(ifftshift(Gc.*conj(Psic)))))/xwidth; % call
VcF = interp1(S,c,S0,'spline');
p = exp(-r*T).*real(fftshift(fft(ifftshift(Gp.*conj(Psip)))))/xwidth; % put
VpF = interp1(S,p,S0,'spline');
cputime_F = toc;
fprintf('%20s%14.10f%14.10f%14.10f\n','Fourier',VcF,VpF,cputime_F)

% Figures
% figures_ft(S,x,xi,Psic,gc,Gc) % call
% figures_ft(S,x,xi,Psip,gp,Gp) % put

%% Monte Carlo

tic;
VcMCb = zeros(nblocks,1);
VpMCb = zeros(nblocks,1);
for i = 1:nblocks

    % Arithmetic Brownian motion X(T) = log(S(T)/S(0)) at time T
    X = muABM*T + sigma*randn(1,nsample)*sqrt(T);

    % Transform to geometric Brownian motion S(T) at time T
    S = S0*exp(X(end));

    % Discounted expected payoff
    VcMCb(i) = exp(-r*T)*mean(max(S-K,0));
    VpMCb(i) = exp(-r*T)*mean(max(K-S,0));

end
VcMC = mean(VcMCb);
VpMC = mean(VpMCb);
scMC = sqrt(var(VcMCb)/nblocks);
spMC = sqrt(var(VpMCb)/nblocks);
cputime_MC = toc;
fprintf('%20s%14.10f%14.10f%14.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)
fprintf('%20s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)
