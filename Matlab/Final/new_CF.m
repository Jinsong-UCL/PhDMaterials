theta = -1;

% Market parameters
S0 = 1; % spot exchange rate
K = 1;
r_0 = [0.05,0.06]; % spot interest rates r_{i0},r_{j0}

% Contract parameters
T = 1; % maturity
%saK = 1; % strike price

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

% Damping parameter
alpha = -2*theta; % Parseval

% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points

% Grids in real and Fourier space
N = ngrid/2;
B = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

% Auxiliary parameters
a_ij_minus = a_i - a_j;
b_ij_minus = b_i - b_j;
a_ij_plus = a_i + a_j;
b_ij_plus = b_i + b_j;
a_ij_rho = rho_v.*a_ij_minus;
b_ij_rho = rho_r.*b_ij_minus;
a_ij_division = a_ij_plus./a_ij_minus;
b_ij_division = b_ij_plus./b_ij_minus;

% Recchioni and Sun, 2016 page 17 eq. (58)
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
% Sun page 14 eq. (52)
c_v = xi.^2-((alpha^2-alpha*a_ij_division')*ones(1,ngrid)+(2*alpha-a_ij_division')*1i*xi);%
% Sun page 15 eq. (53,54)
c_r = xi.^2-((alpha^2-alpha*b_ij_division')*ones(1,ngrid)+(2*alpha-b_ij_division')*1i*xi);%

% Recchioni and Sun, 2016 page 17 eq. (61)
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
% Sun page 16 eq. (60)
d_v = (chi.'*ones(1,ngrid)+(gamma.*a_ij_rho)'.*(1i*xi+alpha));
% Sun page 17 eq. (70,72)
d_r = (lambda.'*ones(1,ngrid)+(eta.*b_ij_rho)'.*(1i*xi+alpha));

% Recchioni and Sun, 2016 page 17 eq. (62)
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
% Sun page 16 eq. (60)
e_v= (d_v.^2 + diag(gamma.^2.*a_ij_minus.^2)*c_v).^0.5;
% Sun page 17 eq. (71,73)
%zeta_qr= 0.5*(4*mu_qr.^2 + 2*diag(eta.^2.*b_ij_minus.^2)*phi_qr).^0.5; %
e_r= (d_r.^2 + diag(eta.^2)*(diag(b_ij_minus.^2)*c_r+2*diag(b_ij_division)*ones(2,1)*(alpha+1i*xi))).^0.5; %

% Recchioni and Sun page, 2016 17 eq. (63)
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
% Sun page 16 eq. (61)
f_v = 1 - exp(-e_v*T);
% Sun page 16 eq. (67)
f_r = 1 - exp(-e_r*T);

% Recchioni and Sun, 2016 page 17 eq. (64)
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
% Sun page 16 eq. (62)
g_v = (e_v-d_v).*exp(-e_v*T)+e_v+d_v;
% Sun page 16 eq. (68)
g_r = (e_r-d_r).*exp(-e_r*T)+e_r+d_r;


% Recchioni and Sun, 2016 page 18 eq. (83)
% Sun page 28 eq. (139)
% W_vq^0
sum_v1 = 2*chi.*v_bar./gamma.^2*log((2*e_v)./g_v);
sum_v2 = chi.*v_bar./gamma.^2*(d_v-e_v)*T;
sum_v3 = v_0./gamma.^2*((d_v.^2-e_v.^2).*f_v./g_v);
phi_v = exp(sum_v1+sum_v2+sum_v3);

% Recchioni and Sun, 2016 page 18 eq. (84)
% Sun page 28 eq. (140)
% W_rq^0
sum_r1 = 2*lambda.*r_bar./eta.^2*log((2*e_r)./g_r);
sum_r2 = lambda.*r_bar./eta.^2*(d_r-e_r)*T;
sum_r3 = r_0./eta.^2*((d_r.^2-e_r.^2).*f_r./g_r);
phi_r = exp(sum_r1+sum_r2+sum_r3);


factor_simple = S0*exp(-r_0(1)*T); % mixes discount and damping
integrand_simple = (K/S0).^(alpha+1+1i*xi).*phi_v.*phi_r./((1i*xi+alpha).*(1i*xi+alpha+1));
priceS_simple = factor_simple*sum(integrand_simple)*dxi/(2*pi);
if theta ==1
    fprintf('The call price of %2.2f is %4.6f\n', K,priceS_simple)
else
    fprintf('The put price of %2.2f is %4.6f\n', K,priceS_simple)
end

% New notation:
fv = xi.^2-((alpha^2-alpha*a_ij_division')*ones(1,ngrid)-(2*alpha-a_ij_division')*1i*xi);%
fr = xi.^2-((alpha^2-alpha*b_ij_division')*ones(1,ngrid)-(2*alpha-b_ij_division')*1i*xi);%
ev = (chi.'*ones(1,ngrid)+(gamma.*a_ij_rho)'.*(-1i*xi+alpha));
er = (lambda.'*ones(1,ngrid)+(eta.*b_ij_rho)'.*(-1i*xi+alpha));
dv= (ev.^2 + diag(gamma.^2.*a_ij_minus.^2)*fv).^0.5;
dr= (er.^2 + diag(eta.^2)*(diag(b_ij_minus.^2)*fr+2*diag(b_ij_division)*ones(2,1)*(alpha-1i*xi))).^0.5; %
gv = (ev - dv)./ (ev + dv);
gr = (er - dr)./ (er + dr);
CFr = exp(lambda.*r_bar./eta.^2*((er-dr)*T-2*log((1-gr.*exp(-dr*T))./(1-gr))) + r_0./eta.^2*((er-dr).*(1-exp(-dr*T))./(1-gr.*exp(-dr*T))));
CFv = exp(chi.*v_bar./gamma.^2*((ev-dv)*T-2*log((1-gv.*exp(-dv*T))./(1-gv))) + v_0./gamma.^2*((ev-dv).*(1-exp(-dv*T))./(1-gv.*exp(-dv*T))));
CF = CFr.*CFv;

factor_simple = S0*exp(-r_0(1)*T); % mixes discount and damping
integrand_new = (K/S0).^(alpha+1-1i*xi).*CF./((-1i*xi+alpha).*(-1i*xi+alpha+1));
priceS_new = factor_simple*sum(integrand_new)*dxi/(2*pi);
if theta ==1
    fprintf('The call price of %2.2f is %4.6f, %4.6f\n', K,priceS_simple, priceS_new)
else
    fprintf('The put price of %2.2f is %4.6f, %4.6f\n', K,priceS_simple, priceS_new)
end
