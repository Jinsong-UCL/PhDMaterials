function [statement] = ECFACF(market,param,fourier,CF_E)
%% Retrieve parameters 
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
xi = fourier.xi; 

% Auxiliary parameters
a_minus = An - Am;
a_plus = An + Am; 

HCF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi(i);
    e1 = kappa + sigma'*rho*a_minus *1i*x;
    a = -Rn + R*1i*x + 0.5*a_plus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    ret = expm([-0.5*e1 -0.5*(sigma'*sigma);a 0.5*e1.']*T);
    B21 = ret(N+1:2*N,1:N);
    B22 = ret(N+1:2*N,N+1:2*N);
    HCF(i) = -2*beta*trace(logm(B22)-0.5*e1.'*T)+trace(B22^(-1)*B21*V_0);
end
HCF_A = exp(HCF);

GCF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi(i);
    e1 = kappa - sigma'*rho*a_minus *1i*x;
    a = -Rn + R*1i*x + 0.5*a_plus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    ret = expm([-0.5*e1 -0.5*(sigma'*sigma);a 0.5*e1.']*T);
    B21 = ret(N+1:2*N,1:N);
    B22 = ret(N+1:2*N,N+1:2*N);
    GCF(i) = -2*beta*trace(logm(B22)-0.5*e1.'*T)+trace(B22^(-1)*B21*V_0);
end
GCF_A = exp(GCF);


figure(1)
plot(xi,real(CF_E),xi,real(HCF_A),xi,real(GCF_A))
axis([-20 20 -0.5 1])
title('Real part of the discounted characteristic function')
xlabel('\xi')
legend('Empirical', 'CF with a plus in E','CF with a minus in E')

figure(2)
plot(xi,imag(CF_E),xi,imag(HCF_A),xi,imag(GCF_A))
axis([-20 20 -0.5 0.5])
title('Imaginary part of the discounted characteristic function')
xlabel('\xi')
legend('Empirical', 'CF with a plus in E','CF with a minus in E')
statement = 1;
end
