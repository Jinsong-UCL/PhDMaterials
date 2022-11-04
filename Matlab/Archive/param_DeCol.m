%% Set parameter values for multifactor stochastic volatility model
% Recchioni, Sun, An explicitly solvable Heston model with stochastic interest rate,
% European Journal of Operational Research 249 (1), 359-377, 2016.

% Market parameters
S0 = 1; % spot exchange rate
r_0 = [0.05,0.06]; % spot interest rates r_{i0},r_{j0}

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Model parameters

param_alpha = 0.5; %
d = 2; % number of volatility factors

% Volatility coefficients or weights
a_i = [0.6650 1.0985];
a_j = [1.6177 1.3588];

% Mean-reversion rate (or strength) of the volatility
chi = [0.9418,1.7909];

% Initial volatility
v_0 = [0.1244,0.0591];

% Long-term average of the volatility
v_bar = [0.037,0.0909];

% Volatility of volatility
gamma = [0.4912,0.08];

% Interest rate coefficients or weights
b_i = [1.004,0.000000];
b_j = [0.000000,1.006];

% Mean-reversion rate (or strength) of the interest rate
lambda = [0.02,0.02]; % lambda_i,lambda_j

% Long-term average of the interest rate
r_bar = [0.05,0.06]; % \bar{r}_i,\bar{r}_j

% Volatility of the interest rate
eta = [0.002,0.002]; % \eta_i,\eta_j

% Correlations
rho_v = [0.5231,-0.398];
rho_r = [-0.23,-0.81];