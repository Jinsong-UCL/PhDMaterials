%% Pricing of European options with the double Heston model and an integral in Fourier space

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algoritm parameters
nsteps = 10; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price
r_i_0 = 0.02; % interest rate

% Mean-reversion rate/strength of the volatility
chi = 0.3;

% Initial volatility
v_0 = 0.05;

% Long-term average of the volatility
v_bar = 0.05;

% Volatility of volatility
gamma = 0.6;

% Relative weight of W^v in SDE of the underlying
delta = 0.0;

% Correlations
rho_v = -0.3;
rho_v_tilde = rho_v + delta;


%%
tic;
% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points
% Grids in real and Fourier space
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

%% Call option
q = 2; % damping factor for a call

% Recchioni Page 17 eq58
% Sun Page 14 eq52
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
phi_vq = 0.5*(xi.^2 + 1i*xi*(2*q-1) - q^2+q);

% Recchioni Page 17 eq61
% Sun Page 17 eq60
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
mu_qv = -0.5*(chi.'*ones(1,ngrid)+(gamma.*rho_v_tilde).'*(1i*xi-q));

% Recchioni Page 17 eq62
% Sun Page 17 eq60
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
zeta_qv= 0.5*(4*mu_qv.^2 + 2*gamma^2*phi_vq*(1+delta^2+2*rho_v*delta)).^0.5;

% Recchioni Page 17 eq63
% Sun Page 17 eq61
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
s_qvg = 1 - exp(-2*zeta_qv*T);

% Recchioni Page 17 eq64
% Sun Page 17 eq62
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
s_qvb = (zeta_qv+mu_qv).*exp(-2*zeta_qv*T)+zeta_qv-mu_qv;

% Recchioni Page 18 eq83
% Sun Page 28 eq139
% W_v_q^0 
sum_1 = 2*chi*v_bar/gamma^2*log(s_qvb/(2*zeta_qv));

sum_2 = 2*chi*v_bar/gamma^2*(mu_qv+zeta_qv)*T;

sum_3 = 2*v_0/gamma^2*(zeta_qv.^2-mu_qv.^2).*s_qvg./s_qvb;

underline_W_q_v = exp(-sum_1-sum_2-sum_3);

% Recchioni Page 6 eq36
% Sun Page 30 eq156
factor = S0*exp(-r_i_0*T*(q-1)); % mixes discount and damping
call_option_integrand = ((S0/K).^(q-1-1i*xi).*exp(-1i*xi*r_i_0*T)).*underline_W_q_v./(-xi.^2-(2*q-1)*xi*1i+q*(q-1));
call_option_price = factor*sum(call_option_integrand)*dxi/(2*pi); 

% 
cputime = toc;
fprintf('%22s%14.10f%14.10f%14.3f\n','Call and put price, CF',call_option_price,0.0,cputime)
% 
% figure(1)
% plot(xi,real(call_option_integrand),xi,imag(call_option_integrand))
% title('Call option integrand')
% xlabel('\xi')
% legend('Real part','Imaginary part')
% figure(2)
% plot(xi,real(put_option_integrand),xi,imag(put_option_integrand))
% title('Put option integrand')
% xlabel('\xi')
% legend('Real part','Imaginary part')