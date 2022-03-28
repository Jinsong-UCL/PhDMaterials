%% Pricing of European options with the double Heston model and an integral in Fourier space
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
a_ij_tilde = a_i + a_j;
rho_v_tilde = rho_v.*a_ij;

%% Call option
q = -2; % damping factor for a call

% Recchioni Page 17 eq58
% Sun Page 14 eq52
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
phi_vq(1,:) = 0.5*(xi.^2 + 1i*xi*(2*q+1) - q^2+q);
phi_vq(2,:) = 0.5*(xi.^2 - 1i*xi*(2*q-1) - q^2+1);

% Recchioni Page 17 eq61
% Sun Page 17 eq60
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
mu_qv = -0.5*(chi.'*ones(1,ngrid)+(gamma.*rho_v_tilde).'*(1i*xi-q));

% Recchioni Page 17 eq62
% Sun Page 17 eq60
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
zeta_qv = zeros(d,ngrid);
zeta_qv(1,:) = 0.5*(4*mu_qv(1,:).^2 + 2*gamma(1)^2*phi_vq(1,:)*a_ij(1)^2).^0.5;
zeta_qv(2,:) = 0.5*(4*mu_qv(2,:).^2 + 2*gamma(2)^2*phi_vq(2,:)*a_ij(2)^2).^0.5;

% Recchioni Page 17 eq63
% Sun Page 17 eq61
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
s_qvg = 1 - exp(-2*zeta_qv*T);

% Recchioni Page 17 eq64
% Sun Page 17 eq62
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
s_qvb = (zeta_qv+mu_qv).*exp(-2*zeta_qv*T)+zeta_qv-mu_qv;

% Recchioni Page 18 eq83
% Sun Page 28 eq139
% W_v_q^0 
sum_1 = 2*chi_1*v_bar_1/gamma_1^2*log(s_qvb(1,:)/(2*zeta_qv(1,:))) ...
        + 2*chi_2*v_bar_2/gamma_2^2*log(s_qvb(2,:)./(2*zeta_qv(2,:)));

sum_2 = 2*chi_1*v_bar_1/gamma_1^2*(mu_qv(1,:)+zeta_qv(1,:))*T ...
        + 2*chi_2*v_bar_2/gamma_2^2*(mu_qv(2,:)+zeta_qv(2,:))*T;

sum_3 = 2*v_1_0/gamma_1^2*(zeta_qv(1,:).^2-mu_qv(1,:).^2).*s_qvg(1,:)./s_qvb(1,:) ...
        + 2*v_2_0/gamma_2^2*(zeta_qv(2,:).^2-mu_qv(2,:).^2).*s_qvg(2,:)./s_qvb(2,:);

underline_W_q_v = exp(-sum_1-sum_2-sum_3);

% Recchioni Page 6 eq36
% Sun Page 30 eq156
factor = S0*exp(-r_i_0*T*(q-1)); % mixes discount and damping
call_option_integrand = ((S0/K).^(q-1-1i*xi).*exp(-1i*xi*r_i_0*T)).*underline_W_q_v./(-xi.^2-(2*q-1)*xi*1i+q*(q-1));
call_option_price = factor*sum(call_option_integrand)*dxi/(2*pi); 



% 
cputime = toc;
fprintf('%22s%14.10f%14.10f%14.3f\n','Call and put price, CF',call_option_price,0.0193775319,cputime)
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