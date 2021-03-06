%% Pricing of European options with the double Heston model and an integral in Fourier space

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Algoritm parameters
nsteps = 10; % monitoring dates
dt = T/nsteps;

% Market parameters
S0 = 1; % spot price

% Mean-reversion rate/strength of the volatility
chi = 0.3;

% Initial volatility
v_0 = 0.05;

% Long-term average of the volatility
v_bar = 0.05;

% Volatility of volatility
gamma = 0.6;

% Initial interest rate
r_0 = 0.02; 

% Mean-reversion rate/strength of the volatility
lambda = 0.3;

% Long-term average of the interest rate
r_bar = 0.05;

% Volatility of interest rate
eta = 0.01;

% Relative weight of W^v in SDE of the underlying
delta = 0.0;

% Relative weight of W^r in SDE of the underlying
omega = 1;

% Correlations
rho_v = -0.3;
rho_r = -0.23;
rho_v_tilde = rho_v + delta;


%%
tic;
% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points
% Grids in real and Fourier space
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

b = xwidth/2; % upper bound of the support in real space
U = S0*exp(b);
L = S0*exp(-b);

% Analytical Fourier transform of the payoff
l = log(L/S0); % lower log barrier
k = log(K/S0); % log strike
u = log(U/S0); % upper log barrier
% Call 
a1 = max(l,k);
b1 = u;
% % Put
% a1 = min(k,u);
% b1 = l;


%% Call option
q = 2; % damping factor for a call

% Recchioni and Sun page 17 eq. (58)
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
phi_q = 0.5*(xi.^2 + 1i*xi*(2*q-1) - q^2+q);

% Recchioni and Sun page 17 eq. (61)
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
mu_qv = -0.5*(chi.'*ones(1,ngrid)+(gamma.*rho_v_tilde).'*(1i*xi-q));

% Recchioni and Sun page 17 eq. (62)
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
zeta_qv= 0.5*(4*mu_qv.^2 + 2*gamma^2*phi_q*(1+delta^2+2*rho_v*delta)).^0.5;

% Recchioni and Sun page 17 eq. (63)
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
s_qvg = 1 - exp(-2*zeta_qv*T);

% Recchioni and Sun page 17 eq. (64)
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
s_qvb = (zeta_qv+mu_qv).*exp(-2*zeta_qv*T)+zeta_qv-mu_qv;


% Recchioni and Sun page 17 eq. (69)
% \mu_{q, r}=-\frac{1}{2}(\lambda+(l k-q) \eta \Omega \rho_{p, r}
mu_qr = -0.5*(lambda.'*ones(1,ngrid)+(eta.*omega*rho_r).'*(1i*xi-q));

% Recchioni and Sun page 17 eq. (70)
% $\zeta_{q,r}=\frac{1}{2}[4\mu_{q,r}^{2}+2 \eta^{2}
% (\varphi_{q}(k)\Omega^{2}-q+\imath k))]^{1/2}$
zeta_qr= 0.5*(4*mu_qr.^2 + 2*eta^2*(phi_q*omega^2-q+1i*xi)).^0.5;

% Recchioni and Sun page 17 eq. (71)
% $s_{q, r, g}=1-e^{-2 \zeta_{q, r} \tau}$
s_qrg = 1 - exp(-2* zeta_qr*T);

% Recchioni and Sun page 17 eq. (72)
% $s_{q,r,b}=(\zeta_{q,r}+\mu_{q,r})e^{-2\zeta_{q,r}\tau}
% +(\zeta_{q, r}-\mu_{q, r})$
s_qrb = (zeta_qr+mu_qr).*exp(-2*zeta_qr*T)+zeta_qr-mu_qr;

% Recchioni and Sun page 18 eq. (83)
% W_vq^0 
sum_v1 = 2*chi*v_bar/gamma^2*log(s_qvb./(2*zeta_qv));
sum_v2 = 2*chi*v_bar/gamma^2*(mu_qv+zeta_qv)*T;
sum_v3 = 2*v_0/gamma^2*(zeta_qv.^2-mu_qv.^2).*s_qvg./s_qvb;
underline_W_vq = exp(-sum_v1-sum_v2-sum_v3);

% Recchioni and Sun page 18 eq. (84)
% W_rq^0 
sum_r1 = 2*lambda*r_bar/eta^2*log(s_qrb./(2*zeta_qr));
sum_r2 = 2*lambda*r_bar/eta^2*(mu_qr+zeta_qr)*T;
sum_r3 = 2*r_0/eta^2*(zeta_qr.^2-mu_qr.^2).*s_qrg./s_qrb;
underline_W_rq = exp(-sum_r1-sum_r2-sum_r3);

% Recchioni and Sun page 18 eq. (80)
% \nu_{r}=2 \lambda \theta / \eta^{2}-1
nu_r = 2*lambda*r_bar/eta^2 - 1;

% M_{q, r}=\frac{2}{\eta^{2}} \frac{s_{q, r, b}}{s_{q, r, g}}
m_qr = 2/eta^2*s_qrb./s_qrg;

% \tilde{r}_{q}=\frac{4 \zeta_{q, r}^{2} r e^{-2 \zeta_{q, r} \tau}}{s_{q, r, b}^{2}}$
tilde_r = 4 * zeta_qr.^2 * r_0 .* exp(-2* zeta_qr *T)./s_qrb.^2;
%% FFT calculation
Psi = exp(-sum_v1-sum_v2-sum_v3-sum_r1-sum_r2-sum_r3);  
xi2 = 1i*xi-q;
G = S0*((exp(b1*(1+xi2))-exp(a1*(1+xi2)))./(1+xi2) ...
    - (exp(k+b1*xi2)-exp(k+a1*xi2))./xi2);
c = exp(-r_0(1)*T).*real(fftshift(fft(ifftshift(G.*conj(Psi)))))/xwidth; % call
VF = interp1(S0*exp(x),c,S0,'spline');

%% Recchioni Calculation
% Recchioni Page 6 eq. (34)
time_factor = T * exp(lambda*T) / (1+ exp(lambda*T));
tail = (m_qr/(m_qr+time_factor))^(nu_r+1).* exp(-time_factor*m_qr.*tilde_r./(m_qr+time_factor));
%tail = 1;
factor = S0*exp(-r_0*T/(1+exp(lambda*T))); % mixes discount and damping
call_integrand = ((S0/K).^(q-1-1i*xi).*exp(-1i*xi*r_0*T)).*underline_W_vq.*underline_W_rq.*tail./(-xi.^2-(2*q-1)*xi*1i+q*(q-1));
call_price = factor*sum(call_integrand)*dxi/(2*pi); 

cputime = toc;
fprintf('%22s%14.7f%14.7f%14.7f%14.3f\n','Call price, MC, Sun, Parseval',0.1177400939,VF,call_price,cputime)
%fprintf('%22s%14.7f%14.7f%14.7f%14.3f\n','Put price, MC, Sun, Parseval',0.0943065201,VF,call_price,cputime)
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