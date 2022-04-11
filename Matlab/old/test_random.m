        
param_alpha = 0; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided
% seed = 100090;
% rng(seed);

% parameters provided by Guido
a_i = [0.2 0.9];
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
nsteps = 50; % number of time steps
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
v_1_0 = 0.05;
v_2_0 = 0.089;
r_i_0 = 0.09; % interest rate
r_j_0 = 0.07; % interest rate

 % Random numbers
rho_v_1 = -0.3;
rho_v_2 = -0.97;
rho_r_i = -0.23;
rho_r_j = -0.81;

% Monte Carlo parameters; 
nblocks = 400;
npaths = 200;
        v_1 = zeros(nsteps,1)
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
            v_2(steps+1) = v_2(steps) + chi_2 * (v_bar_2 - v_2(steps)) * dt + gamma_2 * sqrt(v_2(steps)) * dZ_v_2(steps);
        end
        v_1
        v_2