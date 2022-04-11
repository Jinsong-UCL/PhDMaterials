% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algoritm parameters
nsteps = 10; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
r_i_0 = 0.02; % interest rate

% Mean-reversion rate/strength of the volatility
chi = 0.3;

% Initial volatility
v_0 = 0.05;

% Long-term average of the volatility
v_bar = 0.05;

% Volatility of volatility
gamma = 0.6;

delta = 0.0;

% Correlations
rho_v = -0.3;
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
        v_1 = zeros(nsteps,1);
        v_1 = [v_0;v_1]; 
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
        x = zeros(nsteps+1,1);

        for steps = 1:nsteps
            mu = r_i_0 - 0.5 * v_1(steps) * (1+delta^2+2*rho_v*delta);
            x(steps+1) = x(steps) + dt * mu + v_1(steps)^0.5 * dW_v_1(steps);%
        end
        S_end = S0*exp(x(end));
    
        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);
    
        VcMCb(1,path) = exp(-r_i_0*T)*payoffs_call;
        VpMCb(1,path) = exp(-r_i_0*T)*payoffs_put;
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