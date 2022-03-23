%%%%%%%%%%%%%%%%%%%%%%%%%% Yu Calculation

%% Fourier transform method
tic;
% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points
% Grids in real and Fourier space
tic
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

% operation of parameters 
a_i_j = a_i - a_j;
a_i_j_tilde = a_i + a_j;
rho_v_tilde = rho_v .* a_i_j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation starts here %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Call option %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q =2; % damping factor for a call
% \varphi_{q}^{v_{n}}(k)=\frac{k^{2}}{2}-\frac{1}{2}\left[\left(q^{2}+q \frac{\tilde{a}_{n}^{i, j}}{a_{n}^{i, j}}\right)-\imath k\left(2 q+\frac{\tilde{a}_{n}^{i, j}}{a_{n}^{i, j}}\right)\right]\\
phi_v_q = zeros(d,ngrid);
% phi_v_q(1,1:ngrid) = 0.5*xi.^2 - 0.5 * ((q^2+ q * a_i_j_tilde(1)/a_i_j(1))- 1i * xi *(2*q + a_i_j_tilde(1)/a_i_j(1)));
% phi_v_q(2,1:ngrid) = 0.5*xi.^2 - 0.5 * ((q^2+ q * a_i_j_tilde(2)/a_i_j(2))- 1i * xi *(2*q + a_i_j_tilde(2)/a_i_j(2)));

phi_v_q(1,1:ngrid) = 0.5*xi.^2 - 0.5 * 1i * (2*q-1) * xi - 0.5*(q-1)*q;
phi_v_q(2,1:ngrid) = 0.5*xi.^2 - 0.5 * 1i * (2*q-1) * xi - 0.5*(q-1)*q;

% \mu_{q, v_{n}}=-\frac{1}{2}\left(\chi_{n}+(\imath k-q) \gamma_{n} \tilde{\rho}_{n, v}\right)\\
mu_q_v = zeros(d,ngrid);
mu_q_v(1,1:ngrid) = -0.5.*(chi(1) + (1i* xi -q) * gamma(1) * rho_v_tilde(1));
mu_q_v(2,1:ngrid) = -0.5.*(chi(2) + (1i* xi -q) * gamma(2) * rho_v_tilde(2));

% \zeta_{q, v_{n}}=\frac{1}{2}\left[4 \mu_{q,v_{n}}^{2}+2 \gamma_{n}^{2} \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i, j}\right)^{2}\right]^{1 / 2}\\
zeta_q_v = zeros(d,ngrid);
zeta_q_v(1,1:ngrid) = 0.5 * (4 * mu_q_v(1,1:ngrid).^2 + 2 * gamma(1)^2 * phi_v_q(1,1:ngrid) * a_i_j(1)^2).^0.5;
zeta_q_v(2,1:ngrid) = 0.5 * (4 * mu_q_v(2,1:ngrid).^2 + 2 * gamma(2)^2 * phi_v_q(2,1:ngrid) * a_i_j(2)^2).^0.5;

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
sum_1 = (2 * chi_1 * v_bar_1 / (gamma_1^2)) .* log(s_q_v_b(1,1:ngrid) / (2 * zeta_q_v(1,1:ngrid))) + (2 * chi_2 * v_bar_2 ./ (gamma_2^2)) .* log(s_q_v_b(2,1:ngrid) ./ (2 * zeta_q_v(2,1:ngrid))); 
sum_2 = (2 * chi_1 * v_bar_1 / (gamma_1^2)) .* (mu_q_v(1,1:ngrid)+ zeta_q_v(1,1:ngrid)) .* T + (2 * chi_2 * v_bar_2 / (gamma_2^2)) .* (mu_q_v(2,1:ngrid)+ zeta_q_v(2,1:ngrid)) * T;
sum_3 = (2 * v_1_0 / (gamma_1^2)) * (zeta_q_v(1,1:ngrid).^2 - mu_q_v(1,1:ngrid).^2) .* s_q_v_g(1,1:ngrid) ./ s_q_v_b(1,1:ngrid) + (2 * v_2_0 / (gamma_2^2)) * (zeta_q_v(2,1:ngrid).^2 - mu_q_v(2,1:ngrid).^2) .* s_q_v_g(2,1:ngrid) ./ s_q_v_b(2,1:ngrid);
underline_W_q_v = exp(-sum_1) .* exp(-sum_2) .* exp(-sum_3);


%%%%%%%%%%%%%%%%%%%%%%% call options page 30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call_option_multiplier = exp(-r_i_0 * T) * exp(q * r_i_0 * T) * S0 / (2 * pi) ;
call_option_integrand = ((S0/ K).^(q - 1 - 1i*xi).* exp(-1i*xi*r_i_0*T)) .* underline_W_q_v .* dxi ./ (-xi.^2 - (2*q-1)*xi*1i + q*(q-1));
call_option_integration = sum(call_option_integrand); 
call_option = call_option_multiplier * call_option_integration; 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% put option %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q =-2; % damping factor for a put
% \varphi_{q}^{v_{n}}(k)=\frac{k^{2}}{2}-\frac{1}{2}\left[\left(q^{2}+q \frac{\tilde{a}_{n}^{i, j}}{a_{n}^{i, j}}\right)-\imath k\left(2 q+\frac{\tilde{a}_{n}^{i, j}}{a_{n}^{i, j}}\right)\right]\\
phi_v_q = zeros(d,ngrid);
phi_v_q(1,1:ngrid) = 0.5*xi.^2 - 0.5 * ((q^2+ q * a_i_j_tilde(1)/a_i_j(1))- 1i * xi *(2*q + a_i_j_tilde(1)/a_i_j(1)));
phi_v_q(2,1:ngrid) = 0.5*xi.^2 - 0.5 * ((q^2+ q * a_i_j_tilde(2)/a_i_j(2))- 1i * xi *(2*q + a_i_j_tilde(2)/a_i_j(2)));

% \mu_{q, v_{n}}=-\frac{1}{2}\left(\chi_{n}+(\imath k-q) \gamma_{n} \tilde{\rho}_{n, v}\right)\\
mu_q_v = zeros(d,ngrid);
mu_q_v(1,1:ngrid) = -0.5.*(chi(1) + (1i* xi -q) * gamma(1) * rho_v_tilde(1));
mu_q_v(2,1:ngrid) = -0.5.*(chi(2) + (1i* xi -q) * gamma(2) * rho_v_tilde(2));

% \zeta_{q, v_{n}}=\frac{1}{2}\left[4 \mu_{q,v_{n}}^{2}+2 \gamma_{n}^{2} \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i, j}\right)^{2}\right]^{1 / 2}\\
zeta_q_v = zeros(d,ngrid);
zeta_q_v(1,1:ngrid) = 0.5 * (4 * mu_q_v(1,1:ngrid).^2 + 2 * gamma(1)^2 * phi_v_q(1,1:ngrid) * a_i_j(1)^2).^0.5;
zeta_q_v(2,1:ngrid) = 0.5 * (4 * mu_q_v(2,1:ngrid).^2 + 2 * gamma(2)^2 * phi_v_q(2,1:ngrid) * a_i_j(2)^2).^0.5;

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
sum_1 = (2 * chi_1 * v_bar_1 / (gamma_1^2)) .* log(s_q_v_b(1,1:ngrid) / (2 * zeta_q_v(1,1:ngrid))) + (2 * chi_2 * v_bar_2 ./ (gamma_2^2)) .* log(s_q_v_b(2,1:ngrid) ./ (2 * zeta_q_v(2,1:ngrid))); 
sum_2 = (2 * chi_1 * v_bar_1 / (gamma_1^2)) .* (mu_q_v(1,1:ngrid)+ zeta_q_v(1,1:ngrid)) .* T + (2 * chi_2 * v_bar_2 / (gamma_2^2)) .* (mu_q_v(2,1:ngrid)+ zeta_q_v(2,1:ngrid)) * T;
sum_3 = (2 * v_1_0 / (gamma_1^2)) * (zeta_q_v(1,1:ngrid).^2 - mu_q_v(1,1:ngrid).^2) .* s_q_v_g(1,1:ngrid) ./ s_q_v_b(1,1:ngrid) + (2 * v_2_0 / (gamma_2^2)) * (zeta_q_v(2,1:ngrid).^2 - mu_q_v(2,1:ngrid).^2) .* s_q_v_g(2,1:ngrid) ./ s_q_v_b(2,1:ngrid);
underline_W_q_v = exp(-sum_1) .* exp(-sum_2) .* exp(-sum_3);


%%%%%%%%%%%%%%%%%%%%%%% put options page 30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
put_option_multiplier = exp(-r_i_0 * T) * exp(q * r_i_0 * T) * S0 / (2 * pi) ;
put_option_integrand = ((S0/ K).^(q - 1 - 1i*xi).* exp(-1i*xi*r_i_0*T)) .* underline_W_q_v .* dxi ./ (-xi.^2 - (2*q-1)*xi*1i + q*(q-1));
put_option_integration = sum(put_option_integrand); 
put_option = put_option_multiplier * put_option_integration; 


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