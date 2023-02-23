theta = -1;
% Damping parameter
alpha = -4*theta;

% Market parameters
S0 = 1; % spot exchange rate
r_0 = [0.05,0.06]; % spot interest rates r_{i0},r_{j0}

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Model parameters
param_alpha = 0.5; %
d = 2; % number of volatility factors

% Interest rate coefficients or weights
a_i = [1.004,0.000000];
a_j = [0.000000,1.006];

% Mean-reversion rate (or strength) of the interest rate
kappar = [0.02,0.02]; % lambda_i,lambda_j

% Long-term average of the interest rate
r_bar = [0.05,0.06]; % \bar{r}_i,\bar{r}_j

% Volatility of the interest rate
sigmar = [0.002,0.002]; % \eta_i,\eta_j

% Volatility coefficients or weights
b_i = [0.6650 1.0985];
b_j = [1.6177 1.3588];

% Mean-reversion rate (or strength) of the volatility
kappav = [0.9418,1.7909];
% Initial volatility
v_0 = [0.1244,0.0591];
% Long-term average of the volatility
v_bar = [0.037,0.0909];

% Volatility of volatility
sigmav = [0.4912,0.08];


% Correlations
rho_v = [0.5231,-0.398];
rho_r = [-0.23,-0.81];


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
%xi2 = xi +1i*alpha 

% Auxiliary parameters
a_ij_minus = a_i - a_j;
b_ij_minus = b_i - b_j;
a_ij_plus = a_i + a_j;
b_ij_plus = b_i + b_j;
a_ij_rho = rho_r.*a_ij_minus;
b_ij_rho = rho_v.*b_ij_minus;
a_ij_division = a_ij_plus./a_ij_minus;
b_ij_division = b_ij_plus./b_ij_minus;

% New notation:
fr = xi.^2-((alpha^2-alpha*a_ij_division')*ones(1,ngrid)-(2*alpha-a_ij_division')*1i*xi);%
fv = xi.^2-((alpha^2-alpha*b_ij_division')*ones(1,ngrid)-(2*alpha-b_ij_division')*1i*xi);%
er = (kappar.'*ones(1,ngrid)+(sigmar.*a_ij_rho)'.*(-1i*xi+alpha));
ev = (kappav.'*ones(1,ngrid)+(sigmav.*b_ij_rho)'.*(-1i*xi+alpha));
dr= (er.^2 + diag(sigmar.^2)*(diag(a_ij_minus.^2)*fr+2*diag(a_ij_division)*ones(2,1)*(alpha-1i*xi))).^0.5; %
dv= (ev.^2 + diag(sigmav.^2.*b_ij_minus.^2)*fv).^0.5;
gr = (er - dr)./ (er + dr);
gv = (ev - dv)./ (ev + dv);
CFr = exp(kappar.*r_bar./sigmar.^2*((er-dr)*T-2*log((1-gr.*exp(-dr*T))./(1-gr))) + r_0./sigmar.^2*((er-dr).*(1-exp(-dr*T))./(1-gr.*exp(-dr*T))));
CFv = exp(kappav.*v_bar./sigmav.^2*((ev-dv)*T-2*log((1-gv.*exp(-dv*T))./(1-gv))) + v_0./sigmav.^2*((ev-dv).*(1-exp(-dv*T))./(1-gv.*exp(-dv*T))));
CF = CFr.*CFv;

factor_simple = S0*exp(-r_0(1)*T); % mixes discount and damping
payoff = (K/S0).^(alpha+1+1i*xi)./((1i*xi+alpha).*(1i*xi+alpha+1));
integrand_new = conj(payoff).*CF;
priceS_new = factor_simple*sum(integrand_new)*dxi/(2*pi);
if theta ==1
    fprintf('The call price of %2.2f is %4.6f\n', K, priceS_new)
else
    fprintf('The put price of %2.2f is %4.6f\n', K, priceS_new)
end

%% Barrier options
% Analytical Fourier transform of the payoff
l = -B; % = log(L/C); % lower log barrier
k = log(K/S0); % log strike
u = B; % = log(U/C); % upper log barrier

% Integration bounds
if theta == 1 % call
    a = max(l,k);
    b = u;
else % put
    a = min(k,u);
    b = l;
end

% Green, Fusai, Abrahams 2010 Eq. (3.26) with extension to put option
xi2 = alpha + 1i*xi;
G = S0*((exp(b*(1+xi2))-exp(a*(1+xi2)))./(1+xi2) ...
    - (exp(k+b*xi2)-exp(k+a*xi2))./xi2);
c = S0*exp(-r_0(1)*T).*real(fftshift(fft(ifftshift(G.*conj(CF)))))/xwidth;

priceP = interp1(S0*exp(x),c,S0,'spline')