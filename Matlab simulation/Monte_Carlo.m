%% Compute the YU price of European and lookback options
% The model for the underlying is geometric Brownian motion
% dS = mu*S*dt + sigma*S*dW
% dS/S = (r_0 - r_i)*dt - a_i * diag_v * dW_v - b_i * diag_r^alpha * dW_r
% dv_n = chi_n * (v_n_bar - v_n) * dt + gamma_n * v_n^0.5 * dZ_v_n
% dr_m = lambda_m * (r_m_bar - r_m) * dt + eta_m * r_m^alpha * dZ_r_m

% Model parameter
param_alpha = 0; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided
% seed = 100090;
% rng(seed);

% parameters provided by Guido
a_i = [0.2 0.8];
a_j = [0.3 0.7];
b_i = 0.4;
b_j = 0.6;


% Data below are from Recchioni_2016
r_bar_i = 0.02;
r_bar_j = 0.00044;
v_bar_1 = 0.05;
v_bar_2 = 0.0345;

chi_1 = 0.3;
chi_2 = 0.65;
gamma_1 = 0.6;
gamma_2 = 0.018;
lambda_i = 0.01;
lambda_j = 3.62;
eta_i = 0.01;
eta_j = 0.0098;



% Contract parameters
T = 1; % maturity
K = 1; % strike price
nsteps = 50000; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
v_1_0 = 0.05;
v_2_0 = 0.089;
r_i_0 = 0.09; % interest rate
r_j_0 = 0.09; % interest rate


 % Random numbers
rho_v_1 = -0.3;
rho_v_2 = -0.97;
rho_r_i = -0.23;
rho_r_j = -0.81;


% Monte Carlo parameters; 
npaths = 2000;
MC_result = zeros(nsteps+1,npaths);
Smaxs = zeros(1,npaths);
Smins = zeros(1,npaths);
VcMCb = zeros(1,npaths);
VpMCb = zeros(1,npaths);


%% Monte Carlo

tic;
for i = 1:npaths    
    % Increments of the arithmetic Brownian motion X(t) = log(S(t)/S(0))
    % dX = muABM*dt + sigma*sqrt(dt)*randn(ndates,nsample);
    % dS/S = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r    

    % Generate a random path using the above model
    % Volatility part
    % Model variables
    v_1 = zeros(nsteps,1);
    v_1 = [v_1_0;v_1]; 
    v_2 = zeros(nsteps,1);
    v_2 = [v_2_0;v_2];
    % corr (dW_v_1, dZ_v_1) = rho_v_1
    dW_v_i_1 = randn(nsteps,1);
    dW_v_i_help_1 = randn(nsteps,1);
    dZ_v_i_1 = rho_v_1 * dW_v_i_1 + (1-rho_v_1^2)^0.5 * dW_v_i_help_1;   
    % corr (dW_v_2, dZ_v_2) = rho_v_2
    dW_v_j_1 = randn(nsteps,1);
    dW_v_j_help_1 = randn(nsteps,1);
    dZ_v_j_1 = rho_v_2 * dW_v_j_1 + (1-rho_v_2^2)^0.5 * dW_v_j_help_1; 
    % Convert random numbers into Wiener processes
    dW_v_1 = dW_v_i_1 * dt;
    dW_v_2 = dW_v_j_1 * dt;
    dZ_v_1 = dZ_v_i_1 * dt;
    dZ_v_2 = dZ_v_j_1 * dt;

    for steps = 1:nsteps
        v_1(steps+1,1) = v_1(steps,1) + chi_1 * (v_bar_1 - v_1(steps,1)) * dt + gamma_1 * v_1(steps,1)^0.5 * dZ_v_1(steps,1);
        v_2(steps+1,1) = v_2(steps,1) + chi_2 * (v_bar_2 - v_2(steps,1)) * dt + gamma_2 * v_2(steps,1)^0.5 * dZ_v_2(steps,1);
    end

    % Interest rate part
    % Model variables
    r_i = zeros(nsteps,1);
    r_i = [r_i_0;r_i];
    r_j = zeros(nsteps,1);
    r_j = [r_j_0;r_j];
    % corr (dW_r_i, dZ_r_i) = rho_r_i
    dW_r_i_1 = randn(nsteps,1);
    dW_r_i_help_1 = randn(nsteps,1);
    dZ_r_i_1 = rho_r_i * dW_r_i_1 + (1-rho_r_i^2)^0.5 * dW_r_i_help_1;
    % corr (dW_r_j, dZ_r_j) = rho_r_j
    dW_r_j_1 = randn(nsteps,1);
    dW_r_j_help_1 = randn(nsteps,1);
    dZ_r_j_1 = rho_r_j * dW_r_j_1 + (1-rho_r_j^2)^0.5 * dW_r_j_help_1;
    % Convert random numbers into Wiener processes
    dW_r_i = dW_r_i_1 * dt;
    dW_r_j = dW_r_j_1 * dt;
    dZ_r_i = dZ_r_i_1 * dt;
    dZ_r_j = dZ_r_j_1 * dt;
    
    % dr_m = lambda_m * (r_m_bar - r_m) * dt + eta_m * r_m^alpha * dZ_r_m
    for steps = 1:nsteps        
        r_i(steps+1,1) = r_i(steps,1) + lambda_i * (r_bar_i - r_i(steps,1)) * dt + eta_i * r_i(steps,1)^param_alpha * dZ_r_i(steps,1);
        r_j(steps+1,1) = r_j(steps,1) + lambda_j * (r_bar_j - r_j(steps,1)) * dt + eta_j * r_j(steps,1)^param_alpha * dZ_r_j(steps,1);
    end
    
    sum_v_1 = zeros(nsteps,1);
    sum_v_2 = zeros(nsteps,1);
    sum_r_1 = zeros(nsteps,1);
    sum_r_2 = zeros(nsteps,1);
    muYu = zeros(nsteps,1);
    x_i_j = zeros(nsteps+1,1);
    % Valuation for dx
    for steps = 1:nsteps
        % Block for MuYu
        sum_v_1(steps,1) = (a_i(1)^2 - a_j(1)^2) * v_1(steps,1) + (a_i(2)^2 - a_j(2)^2) * v_2(steps,1);
        sum_r_1(steps,1) = b_i^2 * r_i(steps,1)^(2*param_alpha) - b_j^2 * r_j(steps,1)^(2*param_alpha);
        muYu(steps,1) = r_i(steps,1) - r_j(steps,1) + 0.5 * (sum_v_1(steps) + sum_r_1(steps));
        
        % Block for v
        sum_v_2(steps,1) = (a_i(1) - a_j(1)) * v_1(steps,1)^0.5 * dW_v_1(steps,1) + (a_i(2) - a_j(2)) * v_2(steps,1)^0.5 * dW_v_2(steps,1);
        
        % Block for r
        sum_r_2(steps,1) = b_i * r_i(steps,1)^param_alpha * dW_r_i(steps,1) - b_j * r_j(steps,1)^param_alpha * dW_r_j(steps,1);

        % dx = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r 
        x_i_j(steps+1,1) = x_i_j(steps,1) + muYu(steps,1)*dt + sum_v_2(steps,1) + sum_r_2(steps,1);
    end
    
    % Interest rate integration
    r_diff = r_i - r_j;
    discount_factor = 1;
    for steps = 1:nsteps
        discount_factor = discount_factor * (1 + dt * r_diff(steps,1)); 
    end

    % Calculate the Smax and Smin
    Smax = S0*exp(max(x_i_j,[],1));
    Smin = S0*exp(min(x_i_j,[],1));

    payoffs_call = max(Smax - K,0);
    payoffs_put = max(K - Smin,0);

    VcMCb(1,i) = payoffs_call/discount_factor;
    VpMCb(1,i) = payoffs_put/discount_factor;

    Smaxs(1,i) = Smax;
    Smins(1,i) = Smin;
    MC_result(1:nsteps+1,i) = x_i_j;

end
payoffs_call = max(Smaxs - K,0);
payoffs_put = max(K - Smins,0);


plot(MC_result);
VcMC = mean(VcMCb);
VpMC = mean(VpMCb);
cputime_MC = toc;
fprintf('%20s%14.10f%14.10f%14.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Need to be solved how to intergrate the
%%%%%%%%%%%%%%%%%%%%%%%%%%% interest rate ????



%%%%%%%%%%%%%%%%%%%%%%%%%% Yu Calculation
% Fourier parameters
xwidth = 6; % width of the support in real space
ngrid = 2^8; % number of grid points
alpha = -10; % damping factor for a call










%     % Accumulate the increments
%     X = cumsum(dX,1);
% 
%     % Transform extrema to geometric Brownian motion S(t)
%     Smax = S0*exp(max(X,[],1));
%     Smin = S0*exp(min(X,[],1));
% 
%     % Discounted expected payoff
%     VcMCb(i) = exp(-r*T)*mean(max(Smax-K,0));
%     VpMCb(i) = exp(-r*T)*mean(max(K-Smin,0));
% 
% 
% VcMC = mean(VcMCb);
% VpMC = mean(VpMCb);
% scMC = sqrt(var(VcMCb)/nblocks);
% spMC = sqrt(var(VpMCb)/nblocks);
% cputime_MC = toc;
% fprintf('%20s%14.10f%14.10f%14.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)
% fprintf('%20s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)
% 
% 
% 
% %% Analytical solution
% tic
% muABM = r-q-0.5*sigma^2; % drift coefficient of the arithmetic Brownian motion
% d1 = (log(S0/K)+(r-q+0.5*sigma^2)*T)/(sigma*sqrt(T));
% % d2 = d1 - sigma*sqrt(T);
% d2 = (log(S0/K)+(r-q-0.5*sigma^2)*T)/(sigma*sqrt(T));
% Vca = S0*exp(-q*T)*cdf('Normal',d1,0,1) - K*exp(-r*T)*cdf('Normal',d2,0,1);
% Vpa = K*exp(-r*T)*cdf('Normal',-d2,0,1) - S0*exp(-q*T)*cdf('Normal',-d1,0,1);
% % Put-call parity: Vp = Vc + Kexp(-rT) - S0exp(-q*T)
% cputime_a = toc;
% 
% % Analytical solution provided by Matlab's Financial Toolbox
% tic
% [VcaM,VpaM] = blsprice(S0,K,r,T,sigma,q);
% cputime_aM = toc;
% 
% % Print the results
% fprintf('%20s%14s%14s%14s\n','','call','put','CPU_time/s')
% fprintf('%20s%14.10f%14.10f%14.10f\n','BS analytical',Vca,Vpa,cputime_a)
% fprintf('%20s%14.10f%14.10f%14.10f\n','BS analytical Matlab',VcaM,VpaM,cputime_aM)
% 
% if figures ~= 0
% 
%     % Plot the analytical solution
%     [St,t] = meshgrid(0:.05:2,0:0.025:T);
%     d1 = (log(St/K)+(r-q+0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
%     d2 = (log(St/K)+(r-q-0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
% 
%     close all
%     figure(1)
%     Vc = St.*exp(-q*(T-t)).*cdf('Normal',d1,0,1) - K*exp(-r*(T-t)).*cdf('Normal',d2,0,1);
%     Vc(end,:) = max(St(end,:)-K,0);
%     mesh(St,t,Vc)
%     xlabel('S')
%     ylabel('t')
%     zlabel('V')
%     title('Call')
%     view(-30,24)
%     print('-dpng','bsc.png')
%     
%     figure(2)
%     Vp = K*exp(-r*(T-t)).*cdf('Normal',-d2,0,1) - St.*exp(-q*(T-t)).*cdf('Normal',-d1,0,1);
%     Vp(end,:) = max(K-St(end,:),0);
%     mesh(St,t,Vp)
%     xlabel('S')
%     ylabel('t')
%     zlabel('V')
%     title('Put')
%     view(30,24)
%     print('-dpng','bsp.png')
%     
%     % Plot the analytical solution as a function of the log price
%     k = log(K/S0);
%     [xt,t] = meshgrid(-1:.05:1,0:0.025:T);
%     d1 = (xt-k+(r-q+0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
%     d2 = (xt-k+(r-q-0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
%     
%     figure(3)
%     Vc = S0*(exp(xt-q*(T-t)).*cdf('Normal',d1,0,1) - exp(k-r*(T-t)).*cdf('Normal',d2,0,1));
%     Vc(end,:) = S0*max(exp(xt(end,:))-exp(k),0);
%     mesh(xt,t,Vc)
%     xlabel('x')
%     ylabel('t')
%     zlabel('V')
%     title('Call')
%     view(-30,24)
%     print('-dpng','bscx.png')
%     
%     figure(4)
%     Vp = S0*(exp(k-r*(T-t)).*cdf('Normal',-d2,0,1) - exp(xt-q*(T-t)).*cdf('Normal',-d1,0,1));
%     Vp(end,:) = S0*max(exp(k)-exp(xt(end,:)),0);
%     mesh(xt,t,Vp)
%     xlabel('x')
%     ylabel('t')
%     zlabel('V')
%     title('Put')
%     view(30,24)
%     print('-dpng','bspx.png')
% 
% end
% 
% %% Fourier transform method
% 
% % Grids in real and Fourier space
% tic
% N = ngrid/2;
% b = xwidth/2; % upper bound of the support in real space
% dx = xwidth/ngrid; % grid step in real space
% x = dx*(-N:N-1); % grid in real space
% dxi = pi/b; % Nyquist relation: grid step in Fourier space
% xi = dxi*(-N:N-1); % grid in Fourier space
% 
% % Characteristic function at time T
% xia = xi+1i*alpha; % call
% psi = 1i*muABM*xia-0.5*(sigma*xia).^2; % characteristic exponent
% Psic = exp(psi*T); % characteristic function
% xia = xi-1i*alpha; % put
% psi = 1i*muABM*xia-0.5*(sigma*xia).^2; % characteristic exponent
% Psip = exp(psi*T); % characteristic function
% 
% % These functions provide the characteristic functions of 8 Levy processes
% % param = parameters(1,T,T,r,q); % set the parameters editing parameters.m
% % [x,fc,xi,Psic] = kernel(ngrid,-b,b,param,alpha,0,1); % call
% % [x,fp,xi,Psip] = kernel(ngrid,-b,b,param,-alpha,0,1); % put
% 
% % Fourier transform of the payoff
% b = xwidth/2; % upper bound of the support in real space
% U = S0*exp(b);
% L = S0*exp(-b);
% [~,gc,Gc] = payoff(x,xi,alpha,K,L,U,S0,1); % call
% [S,gp,Gp] = payoff(x,xi,-alpha,K,L,U,S0,0); % put
% 
% % Discounted expected payoff computed with the Plancherel theorem
% c = exp(-r*T).*real(fftshift(fft(ifftshift(Gc.*conj(Psic)))))/xwidth; % call
% VcF = interp1(S,c,S0,'spline');
% p = exp(-r*T).*real(fftshift(fft(ifftshift(Gp.*conj(Psip)))))/xwidth; % put
% VpF = interp1(S,p,S0,'spline');
% cputime_F = toc;
% fprintf('%20s%14.10f%14.10f%14.10f\n','Fourier',VcF,VpF,cputime_F)
% 
% % Figures
% % figures_ft(S,x,xi,Psic,gc,Gc) % call
% % figures_ft(S,x,xi,Psip,gp,Gp) % put
