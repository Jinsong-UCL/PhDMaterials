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
a_i = [0.02 0.09];
a_j = [0.03 0.07];
b_i = 0.4;
b_j = 0.6;

% Data below are from Recchioni_2016
r_bar_i = 0.02;
r_bar_j = 0.00044;
r_bar = [r_bar_i,r_bar_j];

v_bar_1 = 0.05;
v_bar_2 = 0.0345;
v_bar = [v_bar_1,v_bar_2];

chi_1 = 0.3;
chi_2 = 0.65;
chi = [chi_1,chi_2];

gamma_1 = 0.6;
gamma_2 = 0.018;
gamma = [gamma_1,gamma_2];

lambda_i = 0.001;
lambda_j = 0.062;
lambda = [lambda_i,lambda_j];

eta_i = 0.001;
eta_j = 0.0098;
eta = [eta_i,eta_j];



% Contract parameters
T = 1; % maturity
K = 1; % strike price
nsteps = 50; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
v_1_0 = 0.05;
v_2_0 = 0.089;
v_0 = [v_1_0,v_2_0];
r_i_0 = 0.09; % interest rate
r_j_0 = 0.07; % interest rate
r_0 = [r_i_0,r_j_0];


 % Random numbers
rho_v_1 = -0.3;
rho_v_2 = -0.97;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0.23;
rho_r_j = -0.81;
rho_r = [rho_r_i,rho_r_j];


% Monte Carlo parameters; 
nblocks = 40;
npaths = 2000;


%% Monte Carlo

tic;
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
for block = 1:nblocks
    MC_result = zeros(nsteps+1,npaths);
    Smaxs = zeros(1,npaths);
    Smins = zeros(1,npaths);
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    for path = 1:npaths    
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
        dW_v_1_1 = randn(nsteps);
        dW_v_1_help_1 = randn(nsteps);
        dZ_v_1_1 = rho_v_1 * dW_v_1_1 + (1-rho_v_1^2)^0.5 * dW_v_1_help_1;   
        % corr (dW_v_2, dZ_v_2) = rho_v_2
        dW_v_2_1 = randn(nsteps,1);
        dW_v_2_help_1 = randn(nsteps,1);
        dZ_v_2_1 = rho_v_2 * dW_v_2_1 + (1-rho_v_2^2)^0.5 * dW_v_2_help_1; 
        % Convert random numbers into Wiener processes
        dW_v_1 = dW_v_1_1 * sqrt(dt);
        dW_v_2 = dW_v_2_1 * sqrt(dt);
        dZ_v_1 = dZ_v_1_1 * sqrt(dt);
        dZ_v_2 = dZ_v_2_1 * sqrt(dt);
    
        for steps = 1:nsteps
            v_1(steps+1) = max(real(v_1(steps) + chi_1 * (v_bar_1 - v_1(steps)) * dt + gamma_1 * sqrt(v_1(steps)) * dZ_v_1(steps)),0);
            v_2(steps+1) = max(real(v_2(steps) + chi_2 * (v_bar_2 - v_2(steps)) * dt + gamma_2 * sqrt(v_2(steps)) * dZ_v_2(steps)),0);
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
        dW_r_i = dW_r_i_1 * sqrt(dt);
        dW_r_j = dW_r_j_1 * sqrt(dt);
        dZ_r_i = dZ_r_i_1 * sqrt(dt);
        dZ_r_j = dZ_r_j_1 * sqrt(dt);
        
        % dr_m = lambda_m * (r_m_bar - r_m) * dt + eta_m * r_m^alpha * dZ_r_m
        for steps = 1:nsteps        
            r_i(steps+1,1) = max(real(r_i(steps,1) + lambda_i * (r_bar_i - r_i(steps,1)) * dt + eta_i * r_i(steps,1)^param_alpha * dZ_r_i(steps,1)),0);
            r_j(steps+1,1) = max(real(r_j(steps,1) + lambda_j * (r_bar_j - r_j(steps,1)) * dt + eta_j * r_j(steps,1)^param_alpha * dZ_r_j(steps,1)),0);
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
        r_diff_bar = r_bar_i - r_bar_j;
    
        % Calculate the Smax and Smin
    
        S_end = S0*exp(x_i_j(nsteps+1,1));
    
        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);
    
        VcMCb(1,path) = exp(-r_diff_bar*T)*payoffs_call;
        VpMCb(1,path) = exp(-r_diff_bar*T)*payoffs_put;

        MC_result(1:nsteps+1,path) = S0*exp(x_i_j);    
    end
    VcMC(block) = mean(VcMCb);
    VpMC(block) = mean(VpMCb);
end
VcMC_result = mean(VcMC);
VpMC_result = mean(VpMC);
scMC = sqrt(var(VcMC)/nblocks);
spMC = sqrt(var(VpMC)/nblocks);

cputime_MC = toc;

fprintf('%22s%14.10f%14.10f\n','Monte Carlo 1st block',VcMC(1),VpMC(1))
fprintf('%22s%14.10f%14.10f\n','Monte Carlo last block',VcMC(end),VpMC(end))
fprintf('%22s%14.10f%14.10f%14.3f\n','Monte Carlo',VcMC_result,VpMC_result,cputime_MC)
fprintf('%22s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)


%%%%%%%%%%%%%%%%%%%%%%%%%%% Need to be solved how to intergrate the
%%%%%%%%%%%%%%%%%%%%%%%%%%% interest rate ???
%%%%%%%%%%%%%%%%%%%%%%%%%%% Downward sloping???
%%%%%%%%%%%%%%%%%%%%%%%%%%% steady result now




