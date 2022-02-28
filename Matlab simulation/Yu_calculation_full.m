%%%%%%%%%%%%%%%%%%%%%%%%%% Yu Calculation

% Model parameter
param_alpha = 0; % or param_alpha = 0.5
d = 2; % dimention of n, need to be decided
% seed = 100090;
% rng(seed);

% parameters provided by Guido
a_i = [0.02 0.09];
a_j = [0.03 0.07];
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

lambda_i = 0.001;
lambda_j = 0.062;
lambda = [lambda_i,lambda_j];

eta_i = 0.001;
eta_j = 0.0098;
eta = [eta_i,eta_j];



% Contract parameters
T = 1; % maturity
K = 1; % strike price
nsteps = 50; % monitoring dates
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
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points
q = -2; % damping factor for a call
% Grids in real and Fourier space
tic;
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation Starts now %%%%%%%%%%%%%

a_i_j = a_i - a_j;
a_i_j_tilde = a_i + a_j;
rho_v_tilde = rho_v .* a_i_j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of Phi %%%%%%%%%%%%%%%%%%
% \varphi_{q}^{v_{n}}(k)=\frac{k^{2}}{2}-\frac{1}{2}\left[\left(q^{2}+q \frac{\tilde{a}_{n}^{i, j}}{a_{n}^{i, j}}\right)-\imath k\left(2 q+\frac{\tilde{a}_{n}^{i, j}}{a_{n}^{i, j}}\right)\right]\\
phi_v_q = zeros(d,ngrid);
phi_v_q(1,1:ngrid) = 0.5*xi.^2 - 0.5 * ((q^2+ q * a_i_j_tilde(1)/a_i_j(1))- 1i * xi *(2*q + a_i_j_tilde(1)/a_i_j(1)));
phi_v_q(2,1:ngrid) = 0.5*xi.^2 - 0.5 * ((q^2+ q * a_i_j_tilde(2)/a_i_j(2))- 1i * xi *(2*q + a_i_j_tilde(2)/a_i_j(2)));

% \varphi_{q}^{r_{i}}(k)=\frac{k^{2}}{2}-\frac{1}{2}\left[\left(q^{2}+q\right)-\imath k(2 q+1)\right]
phi_ri_q = 0.5 * xi.^2 - 0.5*((q^2+q) - 1i* xi * (2*q+1));
% \varphi_{q}^{r_{j}}(k)=\frac{k^{2}}{2}-\frac{1}{2}\left[\left(q^{2}-q\right)-\imath k(2 q-1)\right]
phi_rj_q = 0.5 * xi.^2 - 0.5*((q^2-q) - 1i* xi * (2*q-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of mu %%%%%%%%%%%%%%%%%%
% \mu_{q, v_{n}}=-\frac{1}{2}\left(\chi_{n}+(\imath k-q) \gamma_{n} \tilde{\rho}_{n, v}\right)\\
mu_q_v = zeros(d,ngrid);
mu_q_v(1,1:ngrid) = -0.5.*(chi(1) + (1i* xi -q) * gamma(1) * rho_v_tilde(1));
mu_q_v(2,1:ngrid) = -0.5.*(chi(2) + (1i* xi -q) * gamma(2) * rho_v_tilde(2));

% \mu_{q, r_{i}}=-\frac{1}{2}\left(\lambda_{i}+(\imath k-q) \eta_{i} \rho_{i, r} b_{i}\right)
mu_q_ri = -0.5 * (lambda_i + (1i * xi - q) * eta_i * rho_r_i * b_i);

% \mu_{q, r_{j}}=-\frac{1}{2}\left(\lambda_{j}+(q-\imath k) \eta_{i} \rho_{i, r} b_{i}\right)
mu_q_rj = -0.5 * (lambda_j + (q - 1i * xi) * eta_j * rho_r_j * b_j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of zeta %%%%%%%%%%%%%%%%%%
% \zeta_{q, v_{n}}=\frac{1}{2}\left[4 \mu_{q,v_{n}}^{2}+2 \gamma_{n}^{2} \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i, j}\right)^{2}\right]^{1 / 2}\\
zeta_q_v = zeros(d,ngrid);
zeta_q_v(1,1:ngrid) = 0.5 * (4 * mu_q_v(1,1:ngrid).^2 + 2 * gamma(1)^2 * phi_v_q(1,1:ngrid) * a_i_j(1)^2).^0.5;
zeta_q_v(2,1:ngrid) = 0.5 * (4 * mu_q_v(2,1:ngrid).^2 + 2 * gamma(2)^2 * phi_v_q(2,1:ngrid) * a_i_j(2)^2).^0.5;

% \zeta_{q, r_{i}}=\frac{1}{2}\left[4 \mu_{q, r_{i}}^{2}+2 \eta_{i}^{2}\left(\varphi_{q}^{r_{i}}(k) b_{i}^{2}-q+\imath k\right)\right]^{1 / 2}
zeta_q_ri = 0.5 * (4 * mu_q_ri.^2) + 2 * eta_i^2 * (phi_ri_q * b_i^2 - q + 1i * k).^0.5;

% \zeta_{q, r_{j}}=\frac{1}{2}\left[4 \mu_{q, r_{j}}^{2}+2 \eta_{j}^{2}\left(\varphi_{q}^{r_{j}}(k) b_{j}^{2}+q-\imath k\right)\right]^{1 / 2}
zeta_q_rj = 0.5 * (4 * mu_q_rj.^2) + 2 * eta_j^2 * (phi_rj_q * b_j^2 + q - 1i * k).^0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of s %%%%%%%%%%%%%%%%%%%%%%%%%%
tau = T - 0;
%s_{q, v_{n}, g} &=1-e^{-2 \zeta_{q, v_{n}} \tau} \\
s_q_v_g = zeros(d,ngrid);
s_q_v_g(1,1:ngrid) = 1 - exp(-2 * zeta_q_v(1,1:ngrid) * tau);
s_q_v_g(2,1:ngrid) = 1 - exp(-2 * zeta_q_v(2,1:ngrid) * tau);

%s_{q, v_{n}, b} &=\left(\zeta_{q, v_{n}}+\mu_{q, v_{n}}\right) e^{-2 \zeta_{q, v_{n}} \tau}+\left(\zeta_{q, v_{n}}-\mu_{q, v_{n}}\right)
s_q_v_b = zeros(d,ngrid);
s_q_v_b(1,1:ngrid) = (zeta_q_v(1,1:ngrid)+ mu_q_v(1,1:ngrid)) .* exp(-2 * zeta_q_v(1,1:ngrid) * tau) + zeta_q_v(1,1:ngrid) - mu_q_v(1,1:ngrid);
s_q_v_b(2,1:ngrid) = (zeta_q_v(2,1:ngrid)+ mu_q_v(2,1:ngrid)) .* exp(-2 * zeta_q_v(2,1:ngrid) * tau) + zeta_q_v(2,1:ngrid) - mu_q_v(2,1:ngrid);






%%%%%%%%%%%%%%%%% W_v_q^0 page 27 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_v_1 = (2 * chi_1 * v_1_0 / (gamma_1^2)) * log(s_q_v_b(1,1:ngrid) / (2 * zeta_q_v(1,1:ngrid))) + (2 * chi_2 * v_2_0 ./ (gamma_2^2)) * log(s_q_v_b(2,1:ngrid) ./ (2 * zeta_q_v(2,1:ngrid)));
sum_v_2 = (2 * chi_1 * v_1_0 / (gamma_1^2)) * (mu_q_v(1,1:ngrid)+ zeta_q_v(1,1:ngrid)) * T + (2 * chi_2 * v_2_0 / (gamma_2^2)) * (mu_q_v(2,1:ngrid)+ zeta_q_v(2,1:ngrid)) * T;
sum_v_3 = (2 * v_1_0 / (gamma_1^2)) * (zeta_q_v(1,1:ngrid).^2 - mu_q_v(1,1:ngrid).^2) .* s_q_v_g(1,1:ngrid) ./ s_q_v_b(1,1:ngrid) + (2 * v_2_0 / (gamma_2^2)) * (zeta_q_v(2,1:ngrid).^2 - mu_q_v(2,1:ngrid).^2) .* s_q_v_g(2,1:ngrid) ./ s_q_v_b(2,1:ngrid);
underline_W_q_v = exp(-sum_v_1) .* exp(-sum_v_2) .* exp(-sum_v_3);


%%%%%%%%%%%%%%%%% W_r_q^0 page 27 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_r_1 = (2 * chi_1 * v_1_0 / (gamma_1^2)) * log(s_q_v_b(1,1:ngrid) / (2 * zeta_q_v(1,1:ngrid))) + (2 * chi_2 * v_2_0 ./ (gamma_2^2)) * log(s_q_v_b(2,1:ngrid) ./ (2 * zeta_q_v(2,1:ngrid)));
sum_r_2 = (2 * chi_1 * v_1_0 / (gamma_1^2)) * (mu_q_v(1,1:ngrid)+ zeta_q_v(1,1:ngrid)) * T + (2 * chi_2 * v_2_0 / (gamma_2^2)) * (mu_q_v(2,1:ngrid)+ zeta_q_v(2,1:ngrid)) * T;
sum_r_3 = (2 * v_1_0 / (gamma_1^2)) * (zeta_q_v(1,1:ngrid).^2 - mu_q_v(1,1:ngrid).^2) .* s_q_v_g(1,1:ngrid) ./ s_q_v_b(1,1:ngrid) + (2 * v_2_0 / (gamma_2^2)) * (zeta_q_v(2,1:ngrid).^2 - mu_q_v(2,1:ngrid).^2) .* s_q_v_g(2,1:ngrid) ./ s_q_v_b(2,1:ngrid);
underline_W_q_r = exp(-sum_r_1) .* exp(-sum_r_2) .* exp(-sum_r_3);



