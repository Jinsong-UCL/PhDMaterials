%% Set parameter values for multifactor stochastic volatility model
% Recchioni, Sun, An explicitly solvable Heston model with stochastic interest rate,
% European Journal of Operational Research 249 (1), 359-377, 2016.

% Model parameters
param_alpha = 0.5; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algoritm parameters
nsteps = 10; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
r_0 = [0.02,0.0]; % initial interest rates r_{i0},r_{j0}

% Volatility coefficients or weights
a_i = [1 1];
a_j = [0.0003 0.0007];

% Mean-reversion rate (or strength) of the volatility
chi = [0.3,0.65];

% Initial volatility
v_0 = [0.05,0.0345];

% Long-term average of the volatility
v_bar = [0.05,0.0345];

% Volatility of volatility
gamma = [0.6,0.018];

% Interest rate coefficients or weights
b_i = [0.0000004,0.000000];
b_j = [0.0000006,0.000000];

% Mean-reversion rate (or strength) of the interest rate
lambda = [0.01,0.02]; % lambda_i,lambda_j

% Long-term average of the interest rate
r_bar = [0.02,0.0]; % \bar{r}_i,\bar{r}_j

% Volatility of the interest rate
eta = [0.001,0.002]; % \eta_i,\eta_j

% Correlations
rho_v = [-0.3,-0.97];
rho_r = [-0.23,-0.81];