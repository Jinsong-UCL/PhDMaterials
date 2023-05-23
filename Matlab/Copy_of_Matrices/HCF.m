function [option_price] = HCF(market,param,fourier,K,theta)
%% Retrieve parameters 
S0 = market.S0;
N = market.d;
T = market.T;

beta = param.beta;
An = param.An;
Am = param.Am;
Rn = param.Rn;
Rm = param.Rm;
R = Rn - Rm;
V_0 = param.V_0;

rho = param.rho;
kappa = param.kappa;
sigma = param.sigma;

ngrid = fourier.ngrid; % number of grid points
dxi = fourier.dxi; 
xi = fourier.xi; 

% Auxiliary parameters
alpha = -2*theta; % Damping parameter
xi_shifted = xi +1i*alpha;
a_minus = An - Am;
a_plus = An + Am; 

CF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi_shifted(i);
    E1 = kappa - rho*sigma*a_minus*1i*x;
    Es = 0.5*(E1+E1.');
    a = -Rn + R*1i*x + 0.5*a_plus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    F = sqrtm(Es*Es - 2*sigma'*sigma* a);
    G = (Es - F)/(Es + F);      
    CF(i) = trace(beta*((Es-F)*T-2*logm((eye(N)-G*expm(-F*T))/(eye(N)-G)))) ...
        + trace(V_0*eye(N)/(sigma'*sigma)*((Es-F)*(eye(N)-expm(-F*T))/(eye(N)-G*expm(-F*T))));
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
end
