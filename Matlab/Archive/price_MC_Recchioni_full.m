%% Pricing of European options with the double Heston model and an integral in Fourier space

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algoritm parameters
nsteps = 10; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price

% Mean-reversion rate/strength of the volatility
chi = 0.3;

% Initial volatility
v_0 = 0.05;

% Long-term average of the volatility
v_bar = 0.05;

% Volatility of volatility
gamma = 0.6;

% Initial interest rate
r_0 = 0.02; 

% Mean-reversion rate/strength of the volatility
lambda = 0.3;

% Long-term average of the interest rate
r_bar = 0.05;

% Volatility of interest rate
eta = 0.01;

% Relative weight of W^v in SDE of the underlying
delta = 0.0;

% Relative weight of W^r in SDE of the underlying
omega = 1;

% Correlations
rho_v = -0.3;
rho_r = -0.23;
rho_v_tilde = rho_v + delta;



% Monte Carlo parameters; 
nblocks = 10000;
npaths = 2000;


%% Monte Carlo

tic;
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    for path = 1:npaths
        v_1 = zeros(nsteps+1,1);
        v_1(1) = v_0;
        % corr (dW_v_1, dZ_v_1) = rho_v_1
        dW_v_1_1 = randn(nsteps,1);
        dW_v_1_help_1 = randn(nsteps,1);
        dZ_v_1_1 = rho_v * dW_v_1_1 + (1-rho_v^2)^0.5 * dW_v_1_help_1;   
        % Convert random numbers into Wiener processes
        dW_v_1 = dW_v_1_1 * sqrt(dt);
        dZ_v_1 = dZ_v_1_1 * sqrt(dt); 

        for steps = 1:nsteps
            v_1(steps+1) = max(real(v_1(steps) + chi * (v_bar - v_1(steps)) * dt + gamma * sqrt(v_1(steps)) * dZ_v_1(steps)),0);
        end
        
        r_1 = zeros(nsteps+1,1);
        r_1(1) = r_0; 
        % corr (dW_v_1, dZ_v_1) = rho_v_1
        dW_r_1_1 = randn(nsteps,1);
        dW_r_1_help_1 = randn(nsteps,1);
        dZ_r_1_1 = rho_r * dW_r_1_1 + (1-rho_r^2)^0.5 * dW_r_1_help_1;   
        % Convert random numbers into Wiener processes
        dW_r_1 = dW_r_1_1 * sqrt(dt);
        dZ_r_1 = dZ_r_1_1 * sqrt(dt); 

        for steps = 1:nsteps
            r_1(steps+1) = max(real(r_1(steps) + lambda * (r_bar - r_1(steps)) * dt + eta * sqrt(r_1(steps)) * dZ_r_1(steps)),0);
        end

        x = zeros(nsteps+1,1);
        for steps = 1:nsteps
            mu = r_1(steps) - 0.5 * (v_1(steps) * (1+delta^2+2*rho_v*delta)+omega^2*r_1(steps));
            x(steps+1) = x(steps) + dt * mu + v_1(steps)^0.5 * dW_v_1(steps) + omega * r_1(steps)^0.5*dW_r_1(steps);%
        end
        S_end = S0*exp(x(end));
    
        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);
    
        VcMCb(1,path) = exp(-r_0*T)*payoffs_call;
        VpMCb(1,path) = exp(-r_0*T)*payoffs_put;
    end
    VcMC(block) = mean(VcMCb);
    VpMC(block) = mean(VpMCb);
end
VcMC_result = mean(VcMC);
VpMC_result = mean(VpMC);
scMC = sqrt(var(VcMC)/nblocks);
spMC = sqrt(var(VpMC)/nblocks);

cputime_MC = toc;
fprintf('%22s%14.10f%14.10f%14.3f\n','Monte Carlo',VcMC_result,VpMC_result,cputime_MC)
fprintf('%22s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)