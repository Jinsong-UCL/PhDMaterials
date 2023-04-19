% Simulation based on our notations 
S0 = 1;
K = 1;
d = 2;
% The following numbers are from Gnoatto and Grasselli 2014
am = [0.7764 0.4837;0.4837 0.9639];
an = [0.6679 0.6277;0.6277 0.8520];
hm = [0.2725 0.0804;0.0804 0.4726];
hn = [0.1841 0.0155;0.0155 0.4761];
y_0 = [0.1688 0.1708;0.1708 0.3169];
% am = [0.7764 0.0;0.0 0.9639];
% an = [0.6679 0.0;0.0 0.8520];
% hm = [0.2725 0.0;0.0 0.4726];
% hn = [0.1841 0.0;0.0 0.4700];
% y_0 = [0.1688 0.0;0.0 0.3169];
%y_bar =
%rho = [-0.5417 0.1899;-0.1170 -0.4834];
%kappa = [1.0426,0.6764;0.9880,0.8778]; 
%sigma = [0.4364,0.1914;0.4966,0.7362];

rho = [-0.5417 0.1899;0.1899 -0.4834];
kappa = [1.0426,0.6764;0.6764,0.8778]; 
sigma = [0.4364,0.1914;0.1914,0.7362];

% rho = [-0.5417 0.0;0.0 -0.4834];
% rho test
if sum(eig(eye(d)-rho*rho.')>=0) == d
    fprintf("rho is valid\n");
end

% kappa = [1.0426,0.0;0.0,0.8778]; 
% sigma = [0.4364,0.0;0.0,0.7362];

h = hm - hn;
beta = 3.1442;
T = 1;

% Number of simulations 
nblocks = 100;
npaths = 100;
% Number of steps 
nsteps = 400
dt = T/nsteps;
% Initialize calls
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);

% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points

% Grids in real and Fourier space
N = ngrid/2;
B = xwidth/2; % upper bound of the support in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

tic;

[msgStr,warnId] = lastwarn;
warnStruct = warning('off',warnId);

phi_e_block = zeros(nblocks,ngrid);% accumulate everything inside of it
for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    phi_e_path = zeros(npaths,ngrid);
    for path = 1:npaths
        interet_rate = zeros(nsteps,1);
        y_latest = y_0;
        x_latest = 0; % log(S0/S0)
        for step = 1:nsteps           
            NW = randn(d,d);
            NB = randn(d,d);
            NZ = NW*rho.' + NB*sqrtm(eye(d)-rho*rho.');
            dW = NW * sqrt(dt);
            dZ = NZ * sqrt(dt);
            % positive semi-definitness
            [V, D] = eig(y_latest);
            D = max(D,0);
            y_latest = V*D*inv(V);

            interet_rate(step) = trace(hm * y_latest);
            % Update X
            sum1 = trace((am - an) * y_latest * (am + an)');
            mu = trace(y_latest * h) +  0.5 * (sum1);
            sum2 = trace((am - an) * sqrtm(y_latest) * dW);
            x_latest = x_latest + mu*dt + sum2;

            % Update V
            y_update = y_latest + (beta*(sigma*sigma') -0.5* kappa*y_latest-0.5* y_latest*kappa) * dt ...
                +0.5*(sigma*sqrtm(y_latest)* dZ + dZ'*sqrtm(y_latest)'*sigma');
            % y_latest = max(y_update,0);
            y_latest = y_update;
        end
        S_end = S0*exp(x_latest);

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);


        % Countinuous Compounding
        r_sum = 0;
        for step = 1:nsteps
            r_sum = r_sum + (interet_rate(step)) *dt;
        end
        dfactor = exp(- r_sum);

        VcMCb(1,path) = dfactor*payoffs_call;
        VpMCb(1,path) = dfactor*payoffs_put;
        phi_e_path(path,:) = exp(1i*xi*x_latest);
    end
    VcMC(block) = mean(VcMCb);
    VpMC(block) = mean(VpMCb);
    phi_e_block(block,:) = mean(phi_e_path);
end
simulated_call = mean(VcMC);
simulated_put = mean(VpMC);
scMC = sqrt(var(VcMC)/nblocks);
spMC = sqrt(var(VpMC)/nblocks);
phi_e = mean(phi_e_block);
phi_e_s = sqrt(var(phi_e_block)/nblocks);

cputime_MC = toc
fprintf('%20s%14.10f%14.10f\n','Monte Carlo',simulated_call,simulated_put)
fprintf('%20s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)

figure(1)
plot(xi,real(phi_e))
axis([-20 20 -0.5 1])
title('Real part of the empirical characteristic function')
xlabel('\xi')
legend('Empirical')

figure(2)
plot(xi,imag(phi_e))
axis([-20 20 -1 1])
title('Imaginary part of the empirical characteristic function')
xlabel('\xi')
legend('Empirical')