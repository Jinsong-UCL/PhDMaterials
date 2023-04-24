function [simulated_call, simulated_put, CF_e] = GGsimulation(market,param)
%% Retrieve parameters 
S0 = market.S0;
K = market.K;
d = market.d;
T = market.T;

beta = param.beta;
An = param.An;
Am = param.Am;
Rn = param.Rn;
Rm = param.Rm;
V_0 = param.V_0;

rho = param.rho;
kappa = param.kappa;
sigma = param.sigma;

R = Rn - Rm;

%% Simulation parameters
% Number of simulations
nblocks = 300;
npaths =  300;
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
CF_e_block = zeros(nblocks,ngrid);% accumulate everything inside of it
for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    phi_e_path = zeros(npaths,ngrid);
    for path = 1:npaths
        interet_rate = zeros(nsteps,1);
        V_latest = V_0;
        x_latest = 0; % log(S0/S0)
        for step = 1:nsteps           
            NW = randn(d,d);
            NB = randn(d,d);
            NZ = NW*rho.' + NB*sqrtm(eye(d)-rho*rho.');
            dW = NW * sqrt(dt);
            dZ = NZ * sqrt(dt);
            % positive semi-definitness
            if sum(eig(V_latest)>=0)<d
                [P, Q] = eig(V_latest);
                Q = max(Q,0);
                V_latest = P*Q/P;
            end

            interet_rate(step) = max(trace(Rn * V_latest),0);
            % Update X
            %sum1 = trace((An - Am) * V_latest * (An + Am));
            sum1 = trace((An - Am) * (An + Am) * V_latest);
            %sum1 = trace(0.5*(An - Am) * V_latest * (An + Am) + 0.5*(An + Am)* V_latest' * (An - Am));
            %sum1 = trace(0.5*(am - an)* (am + an) * y_latest + 0.5* y_latest'*(am + an)* (am - an));
            mu = trace(V_latest * R) +  0.5 * (sum1);
            sum2 = trace((An - Am) * sqrtm(V_latest) * dW); 
            %sum2 = trace(0.5*(An - Am) * sqrtm(V_latest) * dW + 0.5*dW'*sqrtm(V_latest)* (An - Am)');
            x_latest = x_latest + mu*dt + sum2;

            % Update V
            %y_update = y_latest + (beta*(sigma*sigma') -kappa*y_latest) * dt ...
            %    +sigma*sqrtm(y_latest)* dZ; %definitely not working
            V_latest = V_latest + (beta*(sigma'*sigma) -0.5* kappa*V_latest-0.5* V_latest*kappa) * dt ...
                +0.5*(sqrtm(V_latest)* dZ*sigma + sigma'*dZ'*sqrtm(V_latest));
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
    CF_e_block(block,:) = mean(phi_e_path);
end
simulated_call = mean(VcMC);
simulated_put = mean(VpMC);
scMC = sqrt(var(VcMC)/nblocks);
spMC = sqrt(var(VpMC)/nblocks);
CF_e = mean(CF_e_block);
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
