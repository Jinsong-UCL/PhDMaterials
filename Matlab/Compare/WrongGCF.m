function [option_price] = WrongGCF(market,param,fourier,K,theta)
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
    e1 = kappa - sigma'*rho*a_minus *1i*x;
    %a = -Rn + R*1i*x + 0.5*a_plus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    %a = -Rn + R*1i*x + 0.5*a_minus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    a = 0.5*(1i*x*1i*x-1i*x)*a_minus*a_minus + (1i*x-1)*Rn - 1i*x*Rm;
    ret = expm(T*[-0.5*e1 -0.5*(sigma'*sigma);a 0.5*e1.']);
    B21 = ret(N+1:2*N,1:N);
    B22 = ret(N+1:2*N,N+1:2*N);
    CF(i) = -2*beta*trace(logm(B22)-0.5*e1.'*T)+trace(B22^(-1)*B21*V_0);
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
