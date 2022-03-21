%% Set parameter values for multifactor stochastic volatility model
% Recchioni, Sun, An explicitly solvable Heston model with stochastic interest rate,
% European Journal of operational research 249.1 (2016): 359-377.

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algoritm parameters
nsteps = 10; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
r_i_0 = 0.02; % interest rate
r_j_0 = 0.0; % interest rate
r_0 = [r_i_0,r_j_0];

% Model parameters
param_alpha = 0.5; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided

% Volatility coefficient
a_i = [1 1];
a_j = [0.0003 0.0007];

% Mean-reversion rate/strength of the volatility
chi_1 = 0.3;
chi_2 = 0.65;
chi = [chi_1,chi_2];

% Initial volatility
v_1_0 = 0.05;
v_2_0 = 0.0345;
v_0 = [v_1_0,v_2_0];

% Long-term average of the volatility
v_bar_1 = 0.05;
v_bar_2 = 0.0345;
v_bar = [v_bar_1,v_bar_2];

% Volatility of volatility
gamma_1 = 0.6;
gamma_2 = 0.018;
gamma = [gamma_1,gamma_2];

% Interest rate coefficient
b_i = 0.0000004;
b_j = 0.0000006;

% Mean-reversion rate/strength of the interest rate 
lambda_i = 0.01;
lambda_j = 0.02;
lambda = [lambda_i,lambda_j];

% Long-term average of the interest rate
r_bar_i = 0.02;
r_bar_j = 0.000;
r_bar = [r_bar_i,r_bar_j];

% Volatility of the interest rate
eta_i = 0.001;
eta_j = 0.002;
eta = [eta_i,eta_j];

% Correlations
rho_v_1 = -0.3;
rho_v_2 = -0.97;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0.23;
rho_r_j = -0.81;
rho_r = [rho_r_i,rho_r_j];