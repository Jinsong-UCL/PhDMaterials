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
    e1 = kappa - sigma'*rho*a_minus *1i*x;
    %a = -Rn + R*1i*x + 0.5*a_plus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    %a = -Rn + R*1i*x + 0.5*a_minus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    a = 0.5*(1i*x*1i*x-1i*x)*a_minus*a_minus + (1i*x-1)*Rn - 1i*x*Rm;
    ret = expm(T*[-0.5*e1 -0.5*(sigma'*sigma);a 0.5*e1.']);
    B21 = ret(N+1:2*N,1:N);
    B22 = ret(N+1:2*N,N+1:2*N);
    HCF(i) = -2*beta*trace(logm(B22)-0.5*e1.'*T)+trace(B22^(-1)*B21*V_0);
end
HCF_A = exp(HCF);

GCF = zeros(1,ngrid);
for i = 1:ngrid
    x = xi(i);
    e1 = kappa - sigma'*rho*a_minus *1i*x;
    %a = -Rn + R*1i*x + 0.5*a_plus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    a = -Rn + R*1i*x + 0.5*a_minus*a_minus*1i*x - 0.5*a_minus*a_minus*x^2;
    ret = expm(T*[-0.5*e1 -0.5*(sigma'*sigma);a 0.5*e1.']);
    B21 = ret(N+1:2*N,1:N);
    B22 = ret(N+1:2*N,N+1:2*N);
    GCF(i) = -2*beta*trace(logm(B22)-0.5*e1.'*T)+trace(B22^(-1)*B21*V_0);
end
GCF_A = exp(GCF);


figure(1)
plot(xi,real(CF_E),xi,real(GCF_A),"--",xi,real(HCF_A),"--")
axis([-20 20 -0.1 0.7])
title('Real part of the discounted characteristic function for the matrix Heston model')
xlabel('\xi')
legend('Empirical','Analytical\_Corrected','Analytical\_Original')
savefig("Real part of the discounted characteristic function for the matrix Heston model(asymmetric).fig")

figure(2)
plot(xi,imag(CF_E),xi,imag(GCF_A),"--",xi,imag(HCF_A),"--")
axis([-20 20 -0.3 0.3])
title('Imaginary part of the discounted characteristic function for the matrix Heston model')
xlabel('\xi')
legend('Empirical','Analytical\_Corrected','Analytical\_Original')
savefig("Imaginary part of the discounted characteristic function for the matrix Heston model(asymmetric).fig")
statement = 1;
end
