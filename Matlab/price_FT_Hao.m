%% Pricing of European options with the double Heston model and an integral in Fourier space
tic

% Call or put parameter
theta = 1; % 1 for call and -1 for put

% Damping parameter
alpha = -theta*2; % Parseval
q = -alpha; % Sun

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

%L = S0*exp(-B); % lower barrier
%U = S0*exp(B); % upper barrier

% Scale
C = S0;
S = C*exp(x);

% Analytical Fourier transform of the payoff
l = -B; % = log(L/C); % lower log barrier
k = log(K/C); % log strike
u = B; % = log(U/C); % upper log barrier

% Integration bounds
if theta == 1 % call
    a = max(l,k);
    b = u;
else % put
    a = min(k,u);
    b = l;
end

% Green, Fusai, Abrahams 2010 Eq. (3.26) with extension to put option
xi2 = alpha + 1i*xi;
G = C*((exp(b*(1+xi2))-exp(a*(1+xi2)))./(1+xi2) ...
    - (exp(k+b*xi2)-exp(k+a*xi2))./xi2);

% Eliminable discontinuities for xi = 0, otherwise 0/0 = NaN
if (alpha == 0)
    G(floor(end/2)+1) = C*(exp(b)-exp(a)-exp(k)*(b-a));
elseif (alpha == -1)
    G(floor(end/2)+1) = C*(b-a+exp(k-b)-exp(k-a));
end

% Auxiliary parameters
a_ij_minus = a_i - a_j;
b_ij_minus = b_i - b_j;
a_ij_plus = a_i + a_j;
b_ij_plus = b_i + b_j;
a_ij_rho = rho_v.*a_ij_minus;
b_ij_rho = rho_r.*b_ij_minus;
a_ij_division = a_ij_plus./a_ij_minus;
b_ij_division = b_ij_plus./b_ij_minus;

% Recchioni and Sun, 2016 page 17 eq. (58)
% \varphi_q^{v_n}(k) = \frac{k^2}{2}-\frac{1}{2}\left[\left(q^2
% + q\frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)-\imath k\left(2q
% + \frac{\tilde{a}_n^{i,j}}{a_n^{i,j}}\right)\right] \\phi_vq = zeros(d,ngrid);
% Sun page 14 eq. (52)
c_v = ones(d,1)*0.5*xi.^2-0.5*((alpha^2-alpha*a_ij_division')*ones(1,ngrid)+(2*alpha-a_ij_division')*1i*xi);%
% Sun page 15 eq. (53,54)
c_r = ones(2,1)*0.5*xi.^2-0.5*((alpha^2-alpha*b_ij_division')*ones(1,ngrid)+(2*alpha-b_ij_division')*1i*xi);%

% Recchioni and Sun, 2016 page 17 eq. (61)
% \mu_{q,v_n} = -\frac{1}{2}(\chi_n+(\imath k-q)\gamma_n\tilde{\rho}_{n,v})
% Sun page 16 eq. (60)
d_v = -0.5*(chi.'*ones(1,ngrid)+(gamma.*a_ij_rho)'.*(1i*xi+alpha));
% Sun page 17 eq. (70,72)
d_r = -0.5*(lambda.'*ones(1,ngrid)+(eta.*b_ij_rho)'.*(1i*xi+alpha));

% Recchioni and Sun, 2016 page 17 eq. (62)
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
% Sun page 16 eq. (60)
e_v= 0.5*(4*d_v.^2 + 2*diag(gamma.^2.*a_ij_minus.^2)*c_v).^0.5; 
% Sun page 17 eq. (71,73)
%zeta_qr= 0.5*(4*mu_qr.^2 + 2*diag(eta.^2.*b_ij_minus.^2)*phi_qr).^0.5; %
e_r= 0.5*(4*d_r.^2 + 2*diag(eta.^2)*(diag(b_ij_minus.^2)*c_r+diag(b_ij_division)*ones(2,1)*(alpha+1i*xi))).^0.5; %

% Recchioni and Sun page, 2016 17 eq. (63)
% s_{q,v_n,g} = 1-e^{-2\zeta_{q,v_n}\tau}
% Sun page 16 eq. (61)
f_v = 1 - exp(-2*e_v*T);
% Sun page 16 eq. (67)
f_r = 1 - exp(-2*e_r*T);

% Recchioni and Sun, 2016 page 17 eq. (64)
% s_{q,v_n,b} = (\zeta_{q,v_n}+\mu_{q,v_n})e^{-2\zeta_{q,v_n}\tau}
% +\zeta_{q,v_n}-\mu_{q,v_n})
% Sun page 16 eq. (62)
g_v = (e_v+d_v).*exp(-2*e_v*T)+e_v-d_v;
% Sun page 16 eq. (68)
g_r = (e_r+d_r).*exp(-2*e_r*T)+e_r-d_r;


% Recchioni and Sun, 2016 page 18 eq. (83)
% Sun page 28 eq. (139)
% W_vq^0 
sum_v1 = 2*chi.*v_bar./gamma.^2*log(g_v./(2*e_v));
sum_v2 = 2*chi.*v_bar./gamma.^2*(d_v+e_v)*T;
sum_v3 = 2*v_0./gamma.^2*((e_v.^2-d_v.^2).*f_v./g_v); 
underline_W_vq = exp(-sum_v1-sum_v2-sum_v3);

% Recchioni and Sun, 2016 page 18 eq. (84)
% Sun page 28 eq. (140)
% W_rq^0 
sum_r1 = 2*lambda.*r_bar./eta.^2*log(g_r./(2*e_r));
sum_r2 = 2*lambda.*r_bar./eta.^2*(d_r+e_r)*T;
sum_r3 = 2*r_0./eta.^2*((e_r.^2-d_r.^2).*f_r./g_r);
underline_W_rq = exp(-sum_r1-sum_r2-sum_r3);


%% Sun calculation
% Recchioni Page 6, 2016 eq. (34)
% Sun page 31 eq. (154,155)
factor_simple = S0*exp(-r_0(1)*T); % mixes discount and damping
integrand_simple = (K/S0).^(alpha+1+1i*xi).*underline_W_vq.*underline_W_rq./((1i*xi+alpha).*(1i*xi+alpha+1));
priceS_simple = factor_simple*sum(integrand_simple)*dxi/(2*pi); 
cputime = toc;
if theta ==1
    fprintf('%22s%14.10f%14.10f\n','Call price, MC, FTS_Simple, FTS, FTP',VcMC_result,priceS_simple)
    %fprintf('%22s%14.10f%14.10f%14.3f\n','Call price, MC, FT_Sun',VcMC_result,priceS,cputime)
else
    fprintf('%22s%14.10f%14.10f\n','Put price, MC, FTS_Simple, FTS, FTP',VpMC_result,priceS_simple)
    %fprintf('%22s%14.10f%14.10f%14.3f\n','Put  price, MC, FT_Sun',VpMC_result,priceS,cputime)
end

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