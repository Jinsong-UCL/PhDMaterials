function [statement] = ECFACF(market,param,phi_empirical)
%% Retrieve parameters 

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


% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points

% Grids in real and Fourier space
N = ngrid/2;
B = xwidth/2; % upper bound of the support in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

% Auxiliary parameters
a_minus = am - an;
a_plus = am + an; 

CF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi(i);
    E = kappa - a_minus * rho *sigma *1i*x;
    F = a_minus*a_minus*x^2 - a_plus*a_minus*1i*x;
    D = sqrtm(E*E.' + sigma*sigma'* (F - 2 * h*1i*x + 2*hm));
    G = (E - D)/(E + D);    
    CF(i) = trace(beta*((E-D)*T-2*logm((eye(d)-G*expm(-D*T))/(eye(d)-G)))) + trace(y_0*eye(d)/(sigma*sigma')*((E-D)*(eye(d)-expm(-D*T))/(eye(d)-G*expm(-D*T))));
end
phi_analytical = exp(CF);

figure(1)
plot(xi,real(phi_empirical),xi,real(phi_analytical))
axis([-20 20 -0.5 1])
title('Real part of the discounted characteristic function')
xlabel('\xi')
legend('Empirical', 'Analytical')

figure(2)
plot(xi,imag(phi_empirical),xi,imag(phi_analytical))
axis([-20 20 -0.5 0.5])
title('Imaginary part of the discounted characteristic function')
xlabel('\xi')
legend('Empirical', 'Analytical')
statement = 1;
end
