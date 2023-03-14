function option_price = barrierPricing(params,L,U,theta,S0,T,K,r_0)
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
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space
xi_shifted = xi +1i*alpha;

l = log(L/S0); % = log(L/S0); % lower log barrier
k = log(K/S0); % log strike
u = log(U/S0); % = log(U/S0); % upper log barrier

% Integration bounds
if theta == 1 % call
    a = max(l,k);
    b = u;
else % put
    a = min(k,u);
    b = l;
end

% Auxiliary parameters
a_ij_minus = a_i - a_j;
a_ij_plus = a_i + a_j;
a_ij_rho = rho.*a_ij_minus;
a_ij_division = a_ij_plus./a_ij_minus;

% New notation:
f =  xi_shifted.^2-a_ij_division'*1i*xi_shifted;
e = (kappa.'*ones(1,ngrid)+(sigma.*a_ij_rho)'.*(-1i*xi_shifted));
d= (e.^2 + diag(sigma.^2)*(diag(a_ij_minus.^2)*f- 2*diag(h)*ones(4,1)*(1i*xi_shifted))).^0.5; %
g = (e - d)./ (e + d);
CF = exp(kappa.*y_bar./sigma.^2*((e-d)*T-2*log((1-g.*exp(-d*T))./(1-g))) + y_0./sigma.^2*((e-d).*(1-exp(-d*T))./(1-g.*exp(-d*T))));


xi2 = alpha + 1i*xi;
G = S0*((exp(b*(1+xi2))-exp(a*(1+xi2)))./(1+xi2) ...
    - (exp(k+b*xi2)-exp(k+a*xi2))./xi2);
c = S0*exp(-r_0(1)*T).*real(fftshift(fft(ifftshift(G.*conj(CF)))))/xwidth;

option_price = interp1(S0*exp(x),c,S0,'spline');

if theta ==1
    fprintf('The call price of %2.2f is %4.6f\n', K, option_price)
else
    fprintf('The put price of %2.2f is %4.6f\n', K, option_price)
end
end