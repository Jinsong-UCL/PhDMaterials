%% Pricing of European options with the double Heston model and an integral in Fourier space
%%
tic;
% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points
% Grids in real and Fourier space
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

% Auxiliary parameters
a_ij = a_i - a_j;
b_ij = b_i - b_j;
a_ij_tilde = a_i + a_j;
b_ij_tilde = b_i + b_j;
a_ij_rho = rho_v.*a_ij;
b_ij_rho = rho_r.*b_ij;
a_ij_prime = a_ij_tilde ./ a_ij;
b_ij_prime = b_ij_tilde ./ b_ij;
%% Call option
q = -2; % damping factor for a call

% Recchioni and Sun page 17 eq. (58)
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
phi_qv = ones(d,1)*0.5*xi.^2-0.5*((q^2+q*a_ij_prime')*ones(1,ngrid)-(2*q+a_ij_prime')*1i*xi);%
phi_qr = ones(2,1)*0.5*xi.^2-0.5*((q^2+q*b_ij_prime')*ones(1,ngrid)-(2*q+b_ij_prime')*1i*xi);%


% Recchioni and Sun page 17 eq. (61)
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
mu_qv = -0.5*(chi.'*ones(1,ngrid)+(gamma.*a_ij_rho).'*(1i*xi-q));
mu_qr = -0.5*(lambda.'*ones(1,ngrid)+(eta.*b_ij_rho).'*(1i*xi-q));

% Recchioni and Sun page 17 eq. (62)
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
zeta_qv= 0.5*(4*mu_qv.^2 + 2*diag(gamma.^2.*a_ij.^2)*phi_qv).^0.5; 
% zeta_qr= 0.5*(4*mu_qr.^2 + 2*diag(eta.^2.*b_ij)*phi_qr).^0.5; %
zeta_qr= 0.5*(4*mu_qr.^2 + 2*diag(eta.^2)*(diag(b_ij.^2)*phi_qr+diag(b_ij_prime)*ones(2,1)*(-q+1i*xi))).^0.5; %

% Recchioni and Sun page 17 eq. (63)
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
s_qvg = 1 - exp(-2 * zeta_qv*T);
s_qrg = 1 - exp(-2 * zeta_qr*T);

% Recchioni and Sun page 17 eq. (64)
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
s_qvb = (zeta_qv+mu_qv).*exp(-2*zeta_qv*T)+zeta_qv-mu_qv;
s_qrb = (zeta_qr+mu_qr).*exp(-2*zeta_qr*T)+zeta_qr-mu_qr;

nu_v = 2*chi.*v_bar./gamma.^2 - 1;
nu_r = 2*lambda.*r_bar./eta.^2 - 1;

m_qv = diag(2./gamma) * (s_qvb./s_qvg);  
m_qr = diag(2./eta) * (s_qrb./s_qrg);

tilde_v = 4 * diag(v_0) * (zeta_qv.^2 .* exp(-2* zeta_qv *T)./s_qvb.^2);
tilde_r = 4 * diag(r_0) * (zeta_qr.^2 .* exp(-2* zeta_qr *T)./s_qrb.^2);

%%%%%% Characteristic function %%%%
%Psi = 




% Recchioni and Sun page 18 eq. (83)
% W_vq^0 
sum_v1 = 2*chi.*v_bar./gamma.^2*log(s_qvb./(2*zeta_qv));
sum_v2 = 2*chi.*v_bar./gamma.^2*(mu_qv+zeta_qv) * T;
sum_v3 = 2*v_0./gamma.^2 * ((zeta_qv.^2-mu_qv.^2).*s_qvg./s_qvb); 
underline_W_vq = exp(-sum_v1-sum_v2-sum_v3);

% Recchioni and Sun page 18 eq. (84)
% W_rq^0 
sum_r1 = 2*lambda.*r_bar./eta.^2*log(s_qrb./(2*zeta_qr));
sum_r2 = 2*lambda.*r_bar./eta.^2*(mu_qr+zeta_qr)*T;
sum_r3 = 2*r_0./eta.^2*((zeta_qr.^2-mu_qr.^2).*s_qrg./s_qrb);
underline_W_rq = exp(-sum_r1-sum_r2-sum_r3);

% Recchioni Page 6 eq. (34)
time_factor = T * exp(lambda(1)*T) / (1+ exp(lambda(1)*T));
factor = S0*exp(-r_0(1)*T/(1+exp(lambda(1)*T))); % mixes discount and damping
tail = (m_qr(1,:)./(m_qr(1,:)+time_factor)).^(nu_r(1)+1) .* exp(-time_factor*m_qr(1,:).*tilde_r(1,:)./(m_qr(1,:)+time_factor));
call_integrand = ((S0/K).^(q-1-1i*xi).*exp(-1i*xi*r_0(1)*T)).*underline_W_vq.*underline_W_rq.*tail./(-xi.^2-(2*q-1)*xi*1i+q*(q-1));
call_price = factor*sum(call_integrand)*dxi/(2*pi); 

% 
cputime = toc;
% fprintf('%22s%14.10f%14.10f%14.3f\n','Call price, CF',call_price,0.0709481041,cputime)
fprintf('%22s%14.10f%14.10f%14.3f\n','Put price, CF',call_price,0.0149844470,cputime)
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