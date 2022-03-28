%%%%%%%%%%%%%%%%%%%%%%%%%% Yu Calculation

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
zeta_q_ri = 0.5 * (4 * mu_q_ri.^2) + 2 * eta_i^2 * (phi_ri_q * b_i^2 - q + 1i * xi).^0.5;

% \zeta_{q, r_{j}}=\frac{1}{2}\left[4 \mu_{q, r_{j}}^{2}+2 \eta_{j}^{2}\left(\varphi_{q}^{r_{j}}(k) b_{j}^{2}+q-\imath k\right)\right]^{1 / 2}
zeta_q_rj = 0.5 * (4 * mu_q_rj.^2) + 2 * eta_j^2 * (phi_rj_q * b_j^2 + q - 1i * xi).^0.5;


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

% s_{q, r_{m}, g}=1-e^{-2 \zeta_{q, r_{m}} \tau}
s_q_ri_g = 1 - exp(-2 * zeta_q_ri * tau);
s_q_rj_g = 1 - exp(-2 * zeta_q_rj * tau);

% s_{q, r_{m}, b}=\left(\zeta_{q, r_{m}}+\mu_{q, r_{m}}\right) e^{-2 \zeta_{q, r_{m}} \tau}+\left(\zeta_{q, r_{m}}-\mu_{q, r_{m}}\right)
s_q_ri_b = (zeta_q_ri + mu_q_ri) .* exp(-2 * zeta_q_ri) + zeta_q_ri - mu_q_ri;
s_q_rj_b = (zeta_q_rj + mu_q_rj) .* exp(-2 * zeta_q_rj) + zeta_q_rj - mu_q_rj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of M %%%%%%%%%%%%%%%%%%%%%%%%%%
% M_{q, r_{m}}=\frac{2}{\eta_{m}^{2}} \frac{s_{q, r_{m}, b}}{s_{q, r_{m}, g}}
m_q_ri = (2 * s_q_ri_b)./ (eta_i^2 * s_q_ri_g);
m_q_rj = (2 * s_q_rj_b)./ (eta_j^2 * s_q_rj_g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of tilde_r %%%%%%%%%%%%%%%%%%%%%%%%%%
% \tilde{r}_{q, m}=\frac{4\left(\zeta_{q, r_{m}}\right)^{2} r_{m} e^{-2 \zeta_{q, r_{m}} \tau}}{s_{q, r_{m}, b}^{2}}
r_q_i_tilde = (4 * zeta_q_ri.^2 * r_bar_i .* exp(-2 * zeta_q_ri * tau))./(s_q_ri_b.^2);
r_q_j_tilde = (4 * zeta_q_rj.^2 * r_bar_j .* exp(-2 * zeta_q_rj * tau))./(s_q_rj_b.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of nu %%%%%%%%%%%%%%%%%%%%%%%%%%
% \nu_{q, r_{m}}=\left(\lambda \theta_{m}^{Q} / \eta_{m}^{2}\right)-1
nu_q_ri = (lambda_i * r_bar_i)/(eta_i^2) -1;
nu_q_rj = (lambda_j * r_bar_j)/(eta_j^2) -1;

%%%%%%%%%%%%%%%%% W_v_q^0 page 27 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_v_1 = (2 * chi_1 * v_bar_1 / (gamma_1^2)) * log(s_q_v_b(1,1:ngrid) / (2 * zeta_q_v(1,1:ngrid))) + (2 * chi_2 * v_bar_2 ./ (gamma_2^2)) * log(s_q_v_b(2,1:ngrid) ./ (2 * zeta_q_v(2,1:ngrid)));
sum_v_2 = (2 * chi_1 * v_bar_1 / (gamma_1^2)) * (mu_q_v(1,1:ngrid)+ zeta_q_v(1,1:ngrid)) * T + (2 * chi_2 * v_bar_2 / (gamma_2^2)) * (mu_q_v(2,1:ngrid)+ zeta_q_v(2,1:ngrid)) * T;
sum_v_3 = (2 * v_1_0 / (gamma_1^2)) * (zeta_q_v(1,1:ngrid).^2 - mu_q_v(1,1:ngrid).^2) .* s_q_v_g(1,1:ngrid) ./ s_q_v_b(1,1:ngrid) + (2 * v_2_0 / (gamma_2^2)) * (zeta_q_v(2,1:ngrid).^2 - mu_q_v(2,1:ngrid).^2) .* s_q_v_g(2,1:ngrid) ./ s_q_v_b(2,1:ngrid);
underline_W_q_v = exp(-sum_v_1) .* exp(-sum_v_2) .* exp(-sum_v_3);


%%%%%%%%%%%%%%%%% W_r_q^0 page 27 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_r_1 = (2 * lambda_i * r_bar_i / (eta_i^2)) * log(s_q_ri_b / (2 * zeta_q_ri)) + (2 * lambda_j * r_bar_j / (eta_j^2)) * log(s_q_rj_b / (2 * zeta_q_rj));
sum_r_2 = (2 * lambda_i * r_bar_i / (eta_i^2)) * (mu_q_ri + zeta_q_ri) * T + (2 * lambda_j * r_bar_j / (eta_j^2)) * (mu_q_rj + zeta_q_rj) * T;
sum_r_3 = (2 * r_i_0 / (eta_i^2)) * (zeta_q_ri.^2 - mu_q_ri.^2) .* s_q_ri_g ./ s_q_ri_b + (2 * r_j_0 / (eta_j^2)) * (zeta_q_rj.^2 - mu_q_rj.^2) .* s_q_rj_g ./ s_q_rj_b;
underline_W_q_r = exp(-sum_r_1) .* exp(-sum_r_2) .* exp(-sum_r_3);


%%%%%%%%%%%%%%%%%%%%%%% options calculation page 30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_multiplier = (T * exp(lambda_i * T))/(1 + exp(lambda_i * T));
end_multiplier = (m_q_ri./(m_q_ri + t_multiplier)).^(nu_q_ri+1) .* exp(-t_multiplier * (m_q_ri .* r_q_i_tilde)/(m_q_ri + t_multiplier));
option_multiplier_out = exp(-r_bar_i * T) * S0 / (2 * pi);

call_option_integrand = ((S0/ K).^(1 - 1i*xi) .* underline_W_q_v .* underline_W_q_r .* end_multiplier) ./ (-xi.^2 -3*xi *1i + 2) .* dxi;
call_option_integration = sum(call_option_integrand);
call_option = option_multiplier_out * call_option_integration;


put_option_integrand = ((S0/ K).^(-3 - 1i*xi) .* underline_W_q_v .* underline_W_q_r .* end_multiplier) ./ (-xi.^2 +5*xi *1i + 6) .* dxi;
put_option_integration = sum(put_option_integrand);
put_option = option_multiplier_out * put_option_integration;

cputime = toc;
fprintf('%22s%14.10f%14.10f%14.3f\n','Call and put price, CF',call_option,put_option,cputime)

figure(1)
plot(xi,real(call_option_integrand),xi,imag(call_option_integrand))
title('Call option integrand')
xlabel('\xi')
legend('Real part','Imaginary part')

figure(2)
plot(xi,real(put_option_integrand),xi,imag(put_option_integrand))
title('Put option integrand')
xlabel('\xi')
legend('Real part','Imaginary part')
