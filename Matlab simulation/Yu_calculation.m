%%%%%%%%%%%%%%%%%%%%%%%%%% Yu Calculation

% Model parameter
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

lambda_i = 0.01;
lambda_j = 3.62;
lambda = [lambda_i,lambda_j];

eta_i = 0.01;
eta_j = 0.0098;
eta = [eta_i,eta_j];



% Contract parameters
T = 1; % maturity
K = 1; % strike price
nsteps = 5000; % monitoring dates
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


%% Fourier transform method
% Fourier parameters
xwidth = 6; % width of the support in real space
ngrid = 2^8; % number of grid points
q = -2; % damping factor for a call
% Grids in real and Fourier space
tic
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation starts here %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_i_j = a_i - a_j;
a_i_j_tilde = a_i + a_j;
rho_v_tilde = rho_v .* a_i_j;

mu_q_v = -0.5.*(chi + )

