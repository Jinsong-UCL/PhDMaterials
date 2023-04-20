function [simulated_call, simulated_put, phi_e] = GGsimulation(market,param)
%% Retrieve parameters 
S0 = market.S0;
K = market.K;
d = market.d;
T = market.T;

beta = param.beta;
am = param.am;
an = param.an;
hm = param.hm;
hn = param.hn;
y_0 = param.y_0;

rho = param.rho;
kappa = param.kappa;
sigma = param.sigma;

h = hm - hn;

%% Simulation parameters
% Number of simulations
nblocks = 100;
npaths = 100;
% Number of steps
nsteps = 10;
dt = T/nsteps;

%% Fourier parameters
xwidth = 20; % width of the support in real space
ngrid = 2^10; % number of grid points

% Grids in real and Fourier space
N = ngrid/2;
B = xwidth/2; % upper bound of the support in real space
dxi = pi/B; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

tic;

%% Simulation procedure
% Initialization
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
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
            %sum1 = trace((am - an) * (am + an)' * y_latest);
            %sum1 = trace(0.5*(am - an) * y_latest * (am + an)' + 0.5*(am + an)* y_latest' * (am - an)');
            %sum1 = trace(0.5*(am - an)* (am + an)' * y_latest + 0.5* y_latest'*(am + an)* (am - an)' );
            mu = trace(y_latest * h) +  0.5 * (sum1);
            sum2 = trace((am - an) * sqrtm(y_latest) * dW); 
            %sum2 = trace(0.5*(am - an) * sqrtm(y_latest) * dW + 0.5*dW'*sqrtm(y_latest)* (am - an)');
            x_latest = x_latest + mu*dt + sum2;

            % Update V
            %y_update = y_latest + (beta*(sigma*sigma') -kappa*y_latest) * dt ...
            %    +sigma*sqrtm(y_latest)* dZ; %definitely not working
            y_update = y_latest + (beta*(sigma*sigma') -0.5* kappa*y_latest-0.5* y_latest*kappa) * dt ...
                +0.5*(sigma*sqrtm(y_latest)* dZ + dZ'*sqrtm(y_latest)'*sigma');
            
            y_latest = y_update;
        end
        S_end = S0*exp(x_latest);

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);


        % discounting
        r_sum = 0;
        for step = 1:nsteps
            r_sum = r_sum + interet_rate(step) ;
        end
        dfactor = exp(- r_sum*dt);

        VcMCb(1,path) = dfactor*payoffs_call;
        VpMCb(1,path) = dfactor*payoffs_put;
        phi_e_path(path,:) = dfactor*exp(1i*xi*x_latest);
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
%phi_e_s = sqrt(var(phi_e_block)/nblocks);

cputime_MC = toc;
fprintf('%20s%14.10f%14.10f\n','Monte Carlo',simulated_call,simulated_put)
fprintf('%20s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)

% figure(1)
% plot(xi,real(phi_e))
% axis([-20 20 -0.5 1])
% title('Real part of the empirical characteristic function')
% xlabel('\xi')
% legend('Empirical')
% 
% figure(2)
% plot(xi,imag(phi_e))
% axis([-20 20 -0.5 0.5])
% title('Imaginary part of the empirical characteristic function')
% xlabel('\xi')
% legend('Empirical')
end
