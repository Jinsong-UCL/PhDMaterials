function option_price = DCF(params,theta,S0,T,K)
a_i = params.a_i;
a_j = params.a_j;
kappa = params.kappa; 
y_bar = params.y_bar; 
sigma = params.sigma; 
y_0 = params.y_0;
rho = params.rho;
hm = params.hm;
hn = params.hn;
h = hm-hn;

% Damping parameter
alpha = -4*theta;

% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points

% Grids in real and Fourier space
N = ngrid/2;
B = xwidth/2; % upper bound of the support in real space
%dx = xwidth/ngrid; % grid step in real space
%x = dx*(-N:N-1); % grid in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space
xi_shifted = xi +1i*alpha;

% Auxiliary parameters
a_minus = a_i - a_j;
a_plus = a_i + a_j;
a_rho = rho.*a_minus;


% New notation:
%f_c =  diag(a_ij_minus.^2)*(xi_shifted.^2-a_ij_division'*1i*xi_shifted);
%f = diag(a_ij_minus.^2)*ones(4,1)*(xi_shifted.^2) - diag(a_ij_plus.*a_ij_minus)*ones(4,1)*1i*xi_shifted; %plus
f = diag(a_minus.^2)*ones(4,1)*(xi_shifted.^2) - diag(a_plus.*a_minus)*ones(4,1)*1i*xi_shifted; %minus
e = (kappa.'*ones(1,ngrid)+(sigma.*a_rho)'.*(-1i*xi_shifted));
d= (e.^2 + diag(sigma.^2)*(f- 2*diag(h)*ones(4,1)*(1i*xi_shifted) +2 * hm')).^0.5; %
g = (e - d)./ (e + d);
CF = exp(kappa.*y_bar./sigma.^2*((e-d)*T-2*log((1-g.*exp(-d*T))./(1-g))) + y_0./sigma.^2*((e-d).*(1-exp(-d*T))./(1-g.*exp(-d*T))));




CF_new = zeros(1,ngrid);
for i = 1:ngrid
    x = xi_shifted(i);
    E = kappa - rho.*sigma.*a_minus*1i*x;
    al = -hm + h*1i*x + 0.5*a_plus.*a_minus*1i*x - 0.5*a_minus.*a_minus*x^2;
    F = sqrt(E.*E - 2*sigma.*sigma.* al);
    %G = (E - F)/(E + F);      
    A_1 = 2*al.*sinh(F*T/2);
    A_2 = F./y_0.*cosh(F*T/2)+E./y_0.*sinh(F*T/2);
    A = A_1./A_2;
    D = log(F./y_0) + (kappa-F)*T/2 - log((F+E)./(2*y_0) + (F-E)./(2*y_0).*exp(-F*T));
    CF_new(i) = sum(-(kappa.*y_bar*T*1i*x.*a_rho)./sigma - A + 2*D.*kappa.*y_bar./(sigma.*sigma));
end
CF_E = exp(CF_new);

plot(xi, CF-CF_E)

factor_simple = S0;
payoff = (K/S0).^(alpha+1+1i*xi)./((1i*xi+alpha).*(1i*xi+alpha+1));
integrand_new = conj(payoff).*CF;
option_price = factor_simple*sum(integrand_new)*dxi/(2*pi);


