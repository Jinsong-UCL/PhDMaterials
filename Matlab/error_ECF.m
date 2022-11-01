%% Compute the Monte Carlo price of European options with the Yu model
% The model for the underlying is geometric Brownian motion
% dS = mu S dt + sigma S dW
% dS/S = (r_i-r_j)dt - (a_i-a_j)^T diag(v) dW_v - b_i diag(r^alpha) dW_r
% dv_n = chi_n (v_n_bar-v_n) dt + gamma_n \sqrt(v_n) dZ_vn
% dr_m = lambda_m (r_m_bar-r_m) dt + eta_m r_m^alpha dZ_rm

% Call or put parameter
theta = 1; % 1 for call and -1 for put

% Damping parameter
alpha = -2*theta; % Parseval


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

% Algorithm parameters
nsteps = 20;
nblocks = 50;
npaths = 200;

% Monte Carlo
tic
dt = T/nsteps;
phi_e_block = zeros(nblocks,ngrid);% accumulate everything inside of it
for block = 1:nblocks
    phi_e_path = zeros(npaths,ngrid);
    for path = 1:npaths    
        % Increments of the arithmetic Brownian motion X(t) = log(S(t)/S(0))
        % dX = muABM*dt + sigma*sqrt(dt)*randn(ndates,nsample);
        % dS/S = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r    
    
        % Generate a random path using the above model
        % Volatility part
        % Model variables
        v = zeros(nsteps+1,d);
        v(1,:) = v_0;
        % corr (dW_v, dZ_v) = rho_v
        dW_v_1 = randn(nsteps,d);
        dW_v_2 = randn(nsteps,d);
        dW_v_3 = dW_v_1*diag(rho_v) + dW_v_2*diag((1-rho_v.^2).^0.5);

        dW_v = dW_v_1 * sqrt(dt);
        dZ_v = dW_v_3 * sqrt(dt);
        
        % dv = chi * (v_bar - v) * dt + gamma * \sqrt(v) * dZ_v
        for steps = 1:nsteps
            v(steps+1,:) = max(v(steps,:) + (v_bar - v(steps,:)) * diag(chi) * dt ...
                + sqrt(v(steps,:)) .* dZ_v(steps,:)* diag(gamma),0);
        end
    
        % Interest rate part
        % Model variables
        r = zeros(nsteps+1,2);
        r(1,:) = r_0;

        % corr (dW_r_i, dZ_r_i) = rho_r_i
        dW_r_1 = randn(nsteps,2);
        dW_r_2 = randn(nsteps,2);
        dW_r_3 = dW_r_1*diag(rho_r) + dW_r_2*diag((1-rho_r.^2).^0.5);

        dW_r = dW_r_1 * sqrt(dt);
        dZ_r = dW_r_3 * sqrt(dt);

        % dr = lambda * (r_bar - r) * dt + eta * r^alpha * dZ_r
        for steps = 1:nsteps   
            r(steps+1,:) = max(r(steps,:) + (r_bar - r(steps,:)) * diag(lambda) * dt ...
                + r(steps,:).^param_alpha .* dZ_r(steps,:)* diag(eta),0);
        end
        
        sum_v_1 = zeros(nsteps,1);
        sum_v_2 = zeros(nsteps,1);
        sum_r_1 = zeros(nsteps,1);
        sum_r_2 = zeros(nsteps,1);
        mu = zeros(nsteps,1);
        x = zeros(nsteps+1,1);
        % Valuation for dx
        for steps = 1:nsteps
            % Block for MuYu
            % sum_v_1(steps) = (a_i(1)^2 - a_j(1)^2) * v_1(steps) + (a_i(2)^2 - a_j(2)^2) * v_2(steps);
            sum_v_1(steps) = v(steps,:) * (a_i.^2' - a_j.^2');
            % sum_r_1(steps) = b_i^2 * r_i(steps)^(2*param_alpha) - b_j^2 * r_j(steps)^(2*param_alpha);
            sum_r_1(steps) = r(steps,:) * (b_i.^2' - b_j.^2');

            mu(steps) = r(steps,1) - r(steps,2) + 0.5 * (sum_v_1(steps) + sum_r_1(steps));
            
            % Block for v
            sum_v_2(steps) = v(steps,:).^0.5 .* dW_v(steps,:) * (a_i'-a_j');
            
            % Block for r
            sum_r_2(steps) = r(steps,:).^param_alpha .* dW_r(steps,:)*(b_i'-b_j');
    
            % dx = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i *
            % diag_r^alpha * dW_r ???????????????
            x(steps+1) = x(steps) + mu(steps)*dt + sum_v_2(steps) + sum_r_2(steps);
        end
    
        % Calculate the emperical characteristic function
    
        phi_e_path(path,:) = exp(1i*xi*x(end)); 
    end
    phi_e_block(block,:) = mean(phi_e_path);
end
phi_e = mean(phi_e_block);
phi_e_s = sqrt(var(phi_e_block)/nblocks);


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
d_v = 0.5*(chi.'*ones(1,ngrid)+(gamma.*a_ij_rho)'.*(1i*xi+alpha));
% Sun page 17 eq. (70,72)
d_r = 0.5*(lambda.'*ones(1,ngrid)+(eta.*b_ij_rho)'.*(1i*xi+alpha));

% Recchioni and Sun, 2016 page 17 eq. (62)
% \zeta_{q,v_{n}}=\frac{1}{2}\left[4\mu_{q,v_{n}}^{2}+2\gamma_{n}^{2}
% \varphi_{q}^{v_{n}}(k)\left(a_{n}^{i,j}\right)^{2}\right]^{1/2} \\
% Sun page 16 eq. (60)
e_v= 0.5*(4*d_v.^2 + 2*diag(gamma.^2.*a_ij_minus.^2)*c_v).^0.5; 
% Sun page 17 eq. (71,73)
e_r= 0.5*(4*d_r.^2 + 2*diag(eta.^2)*(diag(b_ij_minus.^2)*c_r+diag(b_ij_division)*ones(2,1)*(alpha+1i*xi))).^0.5; %
%e_r= 0.5*(4*d_r.^2 + 2*diag(eta.^2.*b_ij_minus.^2)*c_r).^0.5;

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
g_v = (e_v-d_v).*exp(-2*e_v*T)+e_v+d_v;
% Sun page 16 eq. (68)
g_r = (e_r-d_r).*exp(-2*e_r*T)+e_r+d_r;


% Recchioni and Sun, 2016 page 18 eq. (83)
% Sun page 28 eq. (139)
% W_vq^0 
sum_v1 = 2*chi.*v_bar./gamma.^2*log((2*e_v)./g_v);
sum_v2 = 2*chi.*v_bar./gamma.^2*(d_v-e_v)*T;
sum_v3 = 2*v_0./gamma.^2*((d_v.^2-e_v.^2).*f_v./g_v); 
phi_v = exp(sum_v1+sum_v2+sum_v3);

% Recchioni and Sun, 2016 page 18 eq. (84)
% Sun page 28 eq. (140)
% W_rq^0 
sum_r1 = 2*lambda.*r_bar./eta.^2*log((2*e_r)./g_r);
sum_r2 = 2*lambda.*r_bar./eta.^2*(d_r-e_r)*T;
sum_r3 = 2*r_0./eta.^2*((d_r.^2-e_r.^2).*f_r./g_r);
phi_r = exp(sum_r1+sum_r2+sum_r3);

phi_a = phi_r.*phi_v;

MSE = immse(phi_e,phi_a);
RMSE = rmse(phi_e,phi_a);

cputime_MC = toc;


fprintf('%22s%14.10f\n','Mean equare error',MSE)
fprintf('%22s%14.10f\n','Root mean square error',RMSE)


%fprintf('%20s%14.10f%14.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)
%fprintf('%20s%14.10f\n','Monte Carlo stdev',scMC,spMC)

figure(1)
plot(xi,real(phi_a),xi,real(phi_e))
xlim = ([-200 200]);
title('Real part of the analytical and empirical characteristic function')
xlabel('\xi')
legend('Analytical','Empirical')

figure(2)
plot(xi,imag(phi_a),xi,imag(phi_e))
title('Imaginary part of the analytical and empirical characteristic function')
xlabel('\xi')
legend('Analytical','Empirical')

