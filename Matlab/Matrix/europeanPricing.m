function option_price = europeanPricing(params,theta,S0,T,K)

a_i = diag(params.a_i);
a_j = diag(params.a_j);
kappa = diag(params.kappa); 
y_bar = diag(params.y_bar); 
sigma = diag(params.sigma); 
y_0 = diag(params.y_0);
rho = diag(params.rho);
hm = diag(params.hm);
hn = diag(params.hn);
h = hm-hn;

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
a_ij_minus = a_i - a_j;
a_ij_plus = a_i + a_j;
a_ij_rho = rho*a_ij_minus;

% New notation:
E = kappa*ones(4,ngrid)-sigma*a_ij_rho*ones(4,1)*1i*xi_shifted;

F = a_ij_minus*a_ij_minus'*(ones(4,1)*(xi_shifted.^2)) - a_ij_plus*a_ij_minus*ones(4,1)*1i*xi_shifted;

D = sqrt(E.*E + sigma*sigma'* (F - 2 * h*ones(4,1)*1i*xi_shifted + 2*hm*ones(4,ngrid)));

G = (E - D)./ (E + D);

CF_d = kappa*y_bar/(sigma*sigma')*((E-D)*T-2*log((1-G.*exp(-D*T))./(1-G))) + y_0/(sigma*sigma')*((E-D).*(1-exp(-D*T))./(1-G.*exp(-D*T)));
CF = exp(sum(CF_d));


factor_simple = S0;
payoff = (K/S0).^(alpha+1+1i*xi)./((1i*xi+alpha).*(1i*xi+alpha+1));
integrand_new = conj(payoff).*CF;
option_price = factor_simple*sum(integrand_new)*dxi/(2*pi);
if theta ==1
    fprintf('The call price of %2.2f is %4.6f\n', K, option_price)
else
    fprintf('The put price of %2.2f is %4.6f\n', K, option_price)
end
end