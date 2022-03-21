%% Set parameter values for multifactor stochastic volatility model
% Gauthier and Possamai, Efficient simulation of the double Heston model,
% IUP Journal of Computational Mathematics 4 (3), 23-73, 2011

% Model parameters
param_alpha = 0.5; % exponent of the interest rate
d = 2; % number of volatility factors

% Contract parameters
T = 1; % maturity
K = 61.90*0.7; % strike price

% Algorithm parameters
nsteps = 24;
dt = T/nsteps;

% Market parameters
S0 = 61.90; % spot price
r_i_0 = 0.03; % interest rate
r_j_0 = 0.0; % interest rate
r_0 = [r_i_0,r_j_0];

% Volatility coefficient
a_i = [1 1];
a_j = [0.0001 0.0005];

% Mean-reversion rate/strength of the volatility
chi_1 = 0.9;
chi_2 = 1.2;
chi = [chi_1,chi_2];

% Initial volatility
v_1_0 = 0.6;
v_2_0 = 0.7;
v_0 = [v_1_0,v_2_0];

% Long-term average of the volatility
v_bar_1 = 0.1;
v_bar_2 = 0.15;
v_bar = [v_bar_1,v_bar_2];

% Volatility of volatility
gamma_1 = 0.1;
gamma_2 = 0.2;
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
rho_v_1 = -0.5;
rho_v_2 = -0.5;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0;
rho_r_j = -0;
rho_r = [rho_r_i,rho_r_j];