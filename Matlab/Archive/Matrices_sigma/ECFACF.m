function [statement] = ECFACF(market,param,fourier,CF_E)
%% Retrieve parameters 
d = market.d;
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
xi = fourier.xi; 

% Auxiliary parameters
a_minus = An - Am;
a_plus = An + Am; 

HCF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi(i);
    e1 = kappa - sigma*rho*a_minus*1i*x;
    e2 = kappa'- a_minus*rho'*sigma'*1i*x;
    E = 0.5*(e1+e2);
    F = a_minus*a_minus*x^2 - a_plus*a_minus*1i*x;
    D = sqrtm(E*E + 2*sigma*sigma'* (F - 2 * R*1i*x + 2*Rn));
    G = (E - D)/(E + D);    
    HCF(i) = trace(0.5*beta*((E-D)*T-2*logm((eye(d)-G*expm(-D*T))/(eye(d)-G)))) + trace(V_0*eye(d)/(2*sigma*sigma')*((E-D)*(eye(d)-expm(-D*T))/(eye(d)-G*expm(-D*T))));
end
HCF_A = exp(HCF);

GCF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi(i);
    e1 = -kappa + sigma*rho*a_minus *1i*x;
    e2 = kappa' - a_minus*rho'*sigma' *1i*x;
    G = a_minus*a_minus'*(x^2) - a_minus'*a_plus*1i*x - 2 * R*1i*x+2*Rn;
    ret = expm([0.5*e1 -2*(sigma*sigma');-0.5*G 0.5*e2]*T);
    B21 = ret(d+1:2*d,1:d);
    B22 = ret(d+1:2*d,d+1:2*d);
    GCF(i) = -0.5*beta*trace(logm(B22)+0.5*e1*T)+trace(B22^(-1)*B21*V_0);
end
GCF_A = exp(GCF);


figure(1)
plot(xi,real(CF_E),xi,real(HCF_A),xi,real(GCF_A))
axis([-20 20 -0.5 1])
title('Real part of the discounted characteristic function')
xlabel('\xi')
legend('Empirical', 'HAnalytical','GAnalytical')

figure(2)
plot(xi,imag(CF_E),xi,imag(HCF_A),xi,imag(GCF_A))
axis([-20 20 -0.5 0.5])
title('Imaginary part of the discounted characteristic function')
xlabel('\xi')
legend('Empirical', 'HAnalytical','GAnalytical')
statement = 1;
end
