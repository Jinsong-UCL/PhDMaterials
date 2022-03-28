%% Set parameter values for multifactor stochastic volatility model
% Rouah, The Heston model and its extensions in Matlab and C#. 
% John Wiley & Sons, 2013.

% Model parameter
param_alpha = 0.5; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided

% Contract parameters
T = 0.5; % maturity
K = 100; % strike price

% Algorithm parameters
nsteps = 50; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 100; % spot price
r_i_0 = 0.03; % interest rate
r_j_0 = 0.00; % interest rate
r_0 = [r_i_0,r_j_0];

% Volatility coefficient
a_i = [1 1];
a_j = [0.0001 0.0005];

% Mean-reversion rate/strength of the volatility
chi_1 = 5;
chi_2 = 0.9491;
chi = [chi_1,chi_2];

% Initial volatility
v_1_0 = 0.05;
v_2_0 = 0.000;
v_0 = [v_1_0,v_2_0];

% Long-term average of the volatility
v_bar_1 = 0.05;
v_bar_2 = 0.0;
v_bar = [v_bar_1,v_bar_2];

% Volatility of volatility
gamma_1 = 0.0;
gamma_2 = 0.0;
gamma = [gamma_1,gamma_2];

% Interest rate coefficient
b_i = 0.0000;
b_j = 0.0000;

% Mean-reversion rate/strength of the interest rate 
lambda_i = 0; % arbitrary, e.g. 0 or 1
lambda_j = 0; % arbitrary, e.g. 0 or 1
lambda = [lambda_i,lambda_j];

% Long-term average of the interest rate
r_bar_i = r_i_0;
r_bar_j = r_j_0;
r_bar = [r_bar_i,r_bar_j];

% Volatility of the interest rate
eta_i = 0.000; %
eta_j = 0.000; %
eta = [eta_i,eta_j];

% Correlations
rho_v_1 = -0.8;
rho_v_2 = -0;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0;
rho_r_j = -0;
rho_r = [rho_r_i,rho_r_j];