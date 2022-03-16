% Efficient simulation of the double Heston model
% The IUP Journal of Computational Mathematics 4 (3), 23-73

% Model parameter
param_alpha = 0.5; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided
% seed = 100090;
% rng(seed);

% parameters provided by Guido
a_i = [1 1];
a_j = [0 0.000];
b_i = 0.0000;
b_j = 0.0000;


% Data below are from Recchioni_2016
r_bar_i = 0.03;
r_bar_j = 0.0;
r_bar = [r_bar_i,r_bar_j];

% mean reversion rate
chi_1 = 0.9;
chi_2 = 1.2;
chi = [chi_1,chi_2];

% initial volatility
v_bar_1 = 0.1;
v_bar_2 = 0.15;
v_bar = [v_bar_1,v_bar_2];

% volatility of volatility
gamma_1 = 0.1;
gamma_2 = 0.2;
gamma = [gamma_1,gamma_2];

lambda_i = 0.000;
lambda_j = 0.00;
lambda = [lambda_i,lambda_j];

eta_i = 0.000;
eta_j = 0.000;
eta = [eta_i,eta_j];

% Contract parameters
T = 1; % maturity
K = 61.90*0.7; % strike price
nsteps = 24; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 61.90; % spot price
v_1_0 = 0.6;
v_2_0 = 0.7;
v_0 = [v_1_0,v_2_0];
r_i_0 = 0.03; % interest rate
r_j_0 = 0.0; % interest rate
r_0 = [r_i_0,r_j_0];


 % Random numbers
rho_v_1 = -0.5;
rho_v_2 = -0.5;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0;
rho_r_j = -0;
rho_r = [rho_r_i,rho_r_j];