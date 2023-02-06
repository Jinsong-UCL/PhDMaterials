function option_price = barrierPricing(params,L,U,theta,S0,T,K,r_0)

a_i = params.a_i;
a_j = params.a_j;
kappar = params.kappar; 
r_bar = params.r_bar; 
sigmar = params.sigmar; 
b_i = params.b_i;
b_j = params.b_j;
kappav = params.kappav;
v_0 = params.v_0;
v_bar = params.v_bar;
sigmav = params.sigmav;
rho_v = params.rho_v;
rho_r = params.rho_r;

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
b_ij_minus = b_i - b_j;
a_ij_plus = a_i + a_j;
b_ij_plus = b_i + b_j;
a_ij_rho = rho_r.*a_ij_minus;
b_ij_rho = rho_v.*b_ij_minus;
a_ij_division = a_ij_plus./a_ij_minus;
b_ij_division = b_ij_plus./b_ij_minus;

% New notation:
fr = xi_shifted.^2-a_ij_division'*1i*xi_shifted;%
fv = xi_shifted.^2-b_ij_division'*1i*xi_shifted;%
er = (kappar.'*ones(1,ngrid)+(sigmar.*a_ij_rho)'.*(-1i*xi_shifted));
ev = (kappav.'*ones(1,ngrid)+(sigmav.*b_ij_rho)'.*(-1i*xi_shifted));
dr= (er.^2 + diag(sigmar.^2)*(diag(a_ij_minus.^2)*fr+2*diag(a_ij_division)*ones(2,1)*(-1i*xi_shifted))).^0.5; %
dv= (ev.^2 + diag(sigmav.^2.*b_ij_minus.^2)*fv).^0.5;
gr = (er - dr)./ (er + dr);
gv = (ev - dv)./ (ev + dv);
CFr = exp(kappar.*r_bar./sigmar.^2*((er-dr)*T-2*log((1-gr.*exp(-dr*T))./(1-gr))) + r_0./sigmar.^2*((er-dr).*(1-exp(-dr*T))./(1-gr.*exp(-dr*T))));
CFv = exp(kappav.*v_bar./sigmav.^2*((ev-dv)*T-2*log((1-gv.*exp(-dv*T))./(1-gv))) + v_0./sigmav.^2*((ev-dv).*(1-exp(-dv*T))./(1-gv.*exp(-dv*T))));
CF = CFr.*CFv;


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