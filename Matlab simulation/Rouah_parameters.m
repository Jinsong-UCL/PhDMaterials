% Model parameter
param_alpha = 0.5; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided
% seed = 100090;
% rng(seed);

% parameters provided by Guido
a_i = [1 0];
a_j = [0.000 0.000];


% Data below are from Recchioni_2016
r_bar_i = 0.03;
r_bar_j = 0;
r_bar = [r_bar_i,r_bar_j];


chi_1 = 10.7526;
chi_2 = 0.9491;
chi = [chi_1,chi_2];

v_bar_1 = 0.0330;
v_bar_2 = 0.0257;
v_bar = [v_bar_1,v_bar_2];

gamma_1 = 0.3613;
gamma_2 = 0.0517;
gamma = [gamma_1,gamma_2];


% Contract parameters
T = 1; % maturity
K = 100; % strike price
nsteps = 100; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 100; % spot price
v_1_0 = 0.0252;
v_2_0 = 0.0003;
v_0 = [v_1_0,v_2_0];
r_i_0 = 0.03; % interest rate
r_j_0 = 0.00; % interest rate
r_0 = [r_i_0,r_j_0];


 % Random numbers
rho_v_1 = -0.8;
rho_v_2 = -0;
rho_v = [rho_v_1,rho_v_2];
rho_r_i = -0;
rho_r_j = -0;
rho_r = [rho_r_i,rho_r_j];