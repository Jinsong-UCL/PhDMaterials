%% Set parameter values for multifactor stochastic volatility model
% Random set numbers that inspired by Sun
% No Idea how it works

% Model parameter
param_alpha = 0.5; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algorithm parameters
nsteps = 1; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
r_i_0 = 0.09; % interest rate
r_j_0 = 0.07; % interest rate
r_0 = [r_i_0,r_j_0];

% Volatility coefficient
a_i = [1 0.009];
a_j = [0.003 0.007];

% Mean-reversion rate/strength of the volatility
chi_1 = 0.3;
chi_2 = 0.65;
chi = [chi_1,chi_2];

% Initial volatility
v_1_0 = 0.05;
v_2_0 = 0.089;
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
b_i = 0.4;
b_j = 0.6;

% Mean-reversion rate/strength of the interest rate 
lambda_i = 0.0001;
lambda_j = 0.0062;
lambda = [lambda_i,lambda_j];

% Long-term average of the interest rate
r_bar_i = 0.02;
r_bar_j = 0.00044;
r_bar = [r_bar_i,r_bar_j];

% Volatility of the interest rate
eta_i = 0.0001;
eta_j = 0.00098;
eta = [eta_i,eta_j];

% Correlations
rho_v_1 = -0.3;
rho_v_2 = -0.97;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0.23;
rho_r_j = -0.81;
rho_r = [rho_r_i,rho_r_j];