% Adjusted from Vector code
S0 = 1;
K = 1;
theta = -1;
% The following numbers are from Gnoatto and Grasselli 2014
am = [0.7764 0.4837;0.4837 0.9639];
an = [0.6679 0.6277;0.6277 0.8520];
hm = [0.2725 0.0804;0.0804 0.4726];
hn = [0.1841 0.0155;0.0155 0.4761];
y_0 = [0.1688 0.1708;0.1708 0.3169];
% am = [0.7764 0.0;0.0 0.9639];
% an = [0.6679 0.0;0.0 0.8520];
% hm = [0.2725 0.0;0.0 0.4726];
% hn = [0.1841 0.0;0.0 0.4700];
% y_0 = [0.1688 0.0;0.0 0.3169];
%y_bar =
%rho = [-0.5417 0.1899;-0.1170 -0.4834];
%kappa = [1.0426,0.6764;0.9880,0.8778]; 
%sigma = [0.4364,0.1914;0.4966,0.7362];

rho = [-0.5417 0.1899;0.1899 -0.4834];
kappa = [1.0426,0.6764;0.6764,0.8778]; 
sigma = [0.4364,0.1914;0.1914,0.7362];

% rho = [-0.5417 0.0;0.0 -0.4834];
% rho test
if sum(eig(eye(d)-rho*rho.')>=0) == d
    fprintf("rho is valid\n");
end

% kappa = [1.0426,0.0;0.0,0.8778]; 
% sigma = [0.4364,0.0;0.0,0.7362];

h = hm - hn;
beta = 3.1442;
T = 1;
d = 2;

% Damping parameter
alpha = -4*theta;

% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points

% Grids in real and Fourier space
N = ngrid/2;
B = xwidth/2; % upper bound of the support in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space
xi_shifted = xi +1i*alpha;

% Auxiliary parameters
a_ij_minus = am - an;
a_ij_plus = am + an; 
a_ij_rho = rho*a_ij_minus; 

e1 = kappa*ones(d,ngrid) - a_ij_minus * rho *sigma *ones(d,1)*1i*xi_shifted;
e2 = kappa'*ones(d,ngrid) - sigma' *rho' *a_ij_minus' * ones(d,1)*1i*xi_shifted;

E = 0.5*(e1 + e2);
F = a_ij_minus*a_ij_minus'*(ones(d,1)*(xi_shifted.^2)) - a_ij_plus*a_ij_minus'*ones(d,1)*1i*xi_shifted;
D = sqrt(E.*E + sigma*sigma'* (F - 2 * h*ones(d,1)*1i*xi_shifted + 2*hm*ones(d,ngrid)));

G = (E - D)./ (E + D);

CF_d = eye(d)/(sigma*sigma')*((E-D)*T-2*log((1-G.*exp(-D*T))./(1-G))) + y_0/(sigma*sigma')*((E-D).*(1-exp(-D*T))./(1-G.*exp(-D*T)));
%CF = exp(0.25*beta*sum(CF_d));
CF = exp(beta*sum(CF_d));
factor_simple = S0;
payoff = (K/S0).^(alpha+1+1i*xi)./((1i*xi+alpha).*(1i*xi+alpha+1));
integrand_new = conj(payoff).*CF;
option_price = factor_simple*sum(integrand_new)*dxi/(2*pi);
if theta ==1
    fprintf('The call price of K = %2.2f is %4.6f\n', K, option_price)
else
    fprintf('The put price of K = %2.2f is %4.6f\n', K, option_price)
end