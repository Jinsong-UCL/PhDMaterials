function [option_price] = GCF(market,param,theta)
%% Retrieve parameters 
S0 = market.S0;
K = market.K;
d = market.d;
T = market.T;

beta = param.beta;
am = param.am;
an = param.an;
hm = param.hm;
hn = param.hn;
y_0 = param.y_0;

rho = param.rho;
kappa = param.kappa;
sigma = param.sigma;

h = hm - hn;

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

a_minus = am - an;
a_plus = am + an; 

CF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi_shifted(i);
    e1 = -kappa + sigma'*rho*a_minus *1i*x;
    e2 = kappa' - a_minus*rho'*sigma *1i*x;
    G = a_minus*a_minus'*(x^2) - a_minus'*a_plus*1i*x - 2 * h*1i*x+2*hm;
    ret = expm([0.5*e1 -2*(sigma'*sigma);-0.5*G 0.5*e2]*T);
    B21 = ret(d+1:2*d,1:d);
    B22 = ret(d+1:2*d,d+1:2*d);
    CF(i) = -0.5*beta*trace(logm(B22)+0.5*e1*T)+trace(B22^(-1)*B21*y_0);
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


