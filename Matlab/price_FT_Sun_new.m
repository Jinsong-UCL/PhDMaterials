%% Pricing of European options with the double Heston model and an integral in Fourier space
%%
tic;
% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points
% Grids in real and Fourier space
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid; % grid step in real space
x = dx*(-N:N-1); % grid in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

b = xwidth/2; % upper bound of the support in real space
U = S0*exp(b);
L = S0*exp(-b);

% Analytical Fourier transform of the payoff
l = log(L/S0); % lower log barrier
k = log(K/S0); % log strike
u = log(U/S0); % upper log barrier
% Call 
% a1 = max(l,k);
% b1 = u;
% Put
a1 = min(k,u);
b1 = l;


% Auxiliary parameters
a_ij_minus = a_i - a_j;
b_ij_minus = b_i - b_j;
a_ij_plus = a_i + a_j;
b_ij_plus = b_i + b_j;
a_ij_rho = rho_v.*a_ij_minus;
b_ij_rho = rho_r.*b_ij_minus;
a_ij_division = a_ij_plus ./ a_ij_minus;
b_ij_division = b_ij_plus ./ b_ij_minus;

%% Call option
q = -2; % damping factor for a call

% Recchioni and Sun, 2016 page 17 eq. (58)
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
% Sun page 14 eq. (52)
phi_qv = ones(d,1)*0.5*xi.^2-0.5*((q^2+q*a_ij_division')*ones(1,ngrid)-(2*q+a_ij_division')*1i*xi);%
% Sun page 15 eq. (53,54)
phi_qr = ones(2,1)*0.5*xi.^2-0.5*((q^2+q*b_ij_division')*ones(1,ngrid)-(2*q+b_ij_division')*1i*xi);%


% Recchioni and Sun, 2016 page 17 eq. (61)
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
% Sun page 16 eq. (60)
mu_qv = -0.5*(chi.'*ones(1,ngrid)+(gamma.*a_ij_rho).'*(1i*xi-q));
% Sun page 17 eq. (70,72)
mu_qr = -0.5*(lambda.'*ones(1,ngrid)+(eta.*b_ij_rho).'*(1i*xi-q));

% Recchioni and Sun, 2016 page 17 eq. (62)
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
% Sun page 16 eq. (60)
zeta_qv= 0.5*(4*mu_qv.^2 + 2*diag(gamma.^2.*a_ij_minus.^2)*phi_qv).^0.5; 
% Sun page 17 eq. (71,73)
% zeta_qr= 0.5*(4*mu_qr.^2 + 2*diag(eta.^2.*b_ij)*phi_qr).^0.5; %
zeta_qr= 0.5*(4*mu_qr.^2 + 2*diag(eta.^2)*(diag(b_ij_minus.^2)*phi_qr+diag(b_ij_division)*ones(2,1)*(-q+1i*xi))).^0.5; %

% Recchioni and Sun page, 2016 17 eq. (63)
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
% Sun page 16 eq. (61)
s_qvg = 1 - exp(-2 * zeta_qv*T);
% Sun page 16 eq. (67)
s_qrg = 1 - exp(-2 * zeta_qr*T);

% Recchioni and Sun, 2016 page 17 eq. (64)
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
% Sun page 16 eq. (62)
s_qvb = (zeta_qv+mu_qv).*exp(-2*zeta_qv*T)+zeta_qv-mu_qv;
% Sun page 16 eq. (68)
s_qrb = (zeta_qr+mu_qr).*exp(-2*zeta_qr*T)+zeta_qr-mu_qr;

% Sun page 18 eq. (81)
nu_v = 2*chi.*v_bar./gamma.^2 - 1;
% Sun page 19 eq. (83)
nu_r = 2*lambda.*r_bar./eta.^2 - 1;

% Sun page 18 eq. (80)
m_qv = diag(2./gamma) * (s_qvb./s_qvg);  
% Sun page 19 eq. (82)
m_qr = diag(2./eta) * (s_qrb./s_qrg);

% Sun page 18 eq. (80)
tilde_v = 4 * diag(v_0) * (zeta_qv.^2 .* exp(-2* zeta_qv *T)./s_qvb.^2);
% Sun page 19 eq. (82)
tilde_r = 4 * diag(r_0) * (zeta_qr.^2 .* exp(-2* zeta_qr *T)./s_qrb.^2);

% Recchioni and Sun, 2016 page 18 eq. (83)
% Sun page 28 eq. (139)
% W_vq^0 
sum_v1 = 2*chi.*v_bar./gamma.^2*log(s_qvb./(2*zeta_qv));
sum_v2 = 2*chi.*v_bar./gamma.^2*(mu_qv+zeta_qv) * T;
sum_v3 = 2*v_0./gamma.^2 * ((zeta_qv.^2-mu_qv.^2).*s_qvg./s_qvb); 
underline_W_vq = exp(-sum_v1-sum_v2-sum_v3);

% Recchioni and Sun, 2016 page 18 eq. (84)
% Sun page 28 eq. (140)
% W_rq^0 
sum_r1 = 2*lambda.*r_bar./eta.^2*log(s_qrb./(2*zeta_qr));
sum_r2 = 2*lambda.*r_bar./eta.^2*(mu_qr+zeta_qr)*T;
sum_r3 = 2*r_0./eta.^2*((zeta_qr.^2-mu_qr.^2).*s_qrg./s_qrb);
underline_W_rq = exp(-sum_r1-sum_r2-sum_r3);


%% Characteristic function Calculation
% Sun page 19 eq. (85)
Psi = exp(-sum_v1-sum_v2-sum_v3-sum_r1-sum_r2-sum_r3);  
xi2 = 1i*xi-q;
G = S0*((exp(b1*(1+xi2))-exp(a1*(1+xi2)))./(1+xi2) ...
    - (exp(k+b1*xi2)-exp(k+a1*xi2))./xi2);
c = exp(-r_0(1)*T).*real(fftshift(fft(ifftshift(G.*conj(Psi)))))/xwidth; % call
VF = interp1(S0*exp(x),c,S0,'spline');


%% Sun Calculation
% Recchioni Page 6, 2016 eq. (34)
% Sun page 31 eq. (154,155)
time_factor = T * exp(lambda(1)*T) / (1+ exp(lambda(1)*T));
factor = S0*exp(-r_0(1)*T/(1+exp(lambda(1)*T))); % mixes discount and damping
tail = (m_qr(1,:)./(m_qr(1,:)+time_factor)).^(nu_r(1)+1) .* exp(-time_factor*m_qr(1,:).*tilde_r(1,:)./(m_qr(1,:)+time_factor));
call_integrand = ((S0/K).^(q-1-1i*xi).*exp(-1i*xi*r_0(1)*T)).*underline_W_vq.*underline_W_rq.*tail./(-xi.^2-(2*q-1)*xi*1i+q*(q-1));
call_price = factor*sum(call_integrand)*dxi/(2*pi); 

% 
cputime = toc;
%fprintf('%22s%14.10f%14.10f%14.10f%14.3f\n','Call price, MC, FT, FFT',0.0709481041,call_price,VF,cputime)
fprintf('%22s%14.10f%14.10f%14.10f%14.3f\n','Put price, MC, FT, FFT',0.0149844470,call_price,VF,cputime)
% 
% figure(1)
% plot(xi,real(call_option_integrand),xi,imag(call_option_integrand))
% title('Call option integrand')
% xlabel('\xi')
% legend('Real part','Imaginary part')
% figure(2)
% plot(xi,real(put_option_integrand),xi,imag(put_option_integrand))
% title('Put option integrand')
% xlabel('\xi')
% legend('Real part','Imaginary part')