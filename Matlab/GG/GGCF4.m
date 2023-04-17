% The EXACT CF used by Gnoatto and Grasselli
S0 = 1;
K = 1;
theta = 1;
% The following numbers are from Gnoatto and Grasselli 2014
am = [0.7764 0.4837;0.4837 0.9639];
an = [0.6679 0.6277;0.6277 0.8520];
hm = [0.2725 0.0804;0.0804 0.4726];
hn = [0.1841 0.0155;0.0155 0.4761];
h1 = -0.2218;
h2 = -0.1862;
y_0 = [0.1688 0.1708;0.1708 0.3169];
%am = [0.7764 0.0;0.0 0.9639];
%an = [0.6679 0.0;0.0 0.8520];
%hm = [0.2725 0.0;0.0 0.4726];
%hn = [0.1841 0.0;0.0 0.4761];
%y_0 = [0.1688 0.0;0.0 0.3169];
%y_bar =
rho = [-0.5417 0.1899;-0.1170 -0.4834];
kappa1 = [1.0426,0.6764;0.9880,0.8778];
kappa = -0.5*kappa1;
sigma1 = [0.4364,0.1914;0.4966,0.7362];
sigma = 0.5*sigma1;

%rho = [-0.5417 0.1899;0.1899 -0.4834];
%kappa = [1.0426,0.6764;0.6764,0.8778]; 
%sigma = [0.4364,0.1914;0.1914,0.7362];

%rho = [-0.5417 0.0;0.0 -0.4834];
%kappa = [1.0426,0.0;0.0,0.8778]; 
%sigma = [0.4364,0.0;0.0,0.7362];

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
CF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi_shifted(i);
    e1 = kappa + sigma* rho*a_ij_minus *1i*x;
    e2 = kappa' + a_ij_minus'*rho'*sigma' *1i* x;
    D = 0.5*((1i*x)^2-1i*x)*a_ij_minus*(a_ij_minus')+(1i*x-1)*hm-1i*x*hn;
    ret = expm([e1 -2*(sigma'*sigma);D -e2]*T);
    B21 = ret(d+1:2*d,1:d);
    B22 = ret(d+1:2*d,d+1:2*d);
    CF(i) = (1i*x*(h1-h2)-h1)*T -0.5*beta*trace(logm(B22)+e2*T)...
        +trace(B22^(-1)*B21*y_0);
end
CF_E = exp(CF);
factor_simple = S0;
payoff = (K/S0).^(alpha+1+1i*xi)./((1i*xi+alpha).*(1i*xi+alpha+1));
integrand_new = conj(payoff).*CF_E;
option_price = factor_simple*sum(integrand_new)*dxi/(2*pi);
if theta ==1
    fprintf('The call price of %2.2f is %4.6f\n', K, option_price)
else
    fprintf('The put price of %2.2f is %4.6f\n', K, option_price)
end

