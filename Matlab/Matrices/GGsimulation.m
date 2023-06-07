function [simulated_call, simulated_put,scMC,spMC, CF_e] = GGsimulation(market,param,fourier,K,n)
%% Retrieve parameters 
S0 = market.S0;
N = market.d;
T = market.T;

beta = param.beta;
An = param.An;
Am = param.Am;
Rn = param.Rn;
Rm = param.Rm;
hn = param.hn;
hm = param.hm;
R = Rn - Rm;
V_0 = param.V_0;

rho = param.rho;
kappa = param.kappa;
sigma = param.sigma;
sum(eig(sigma)>=0);

ngrid = fourier.ngrid; % number of grid points
xi = fourier.xi; 
%% Simulation parameters
% Number of simulations
nblocks = param.nblocks;
npaths =  param.npaths;
% Number of steps
nsteps = n*50;
dt = T/nsteps;

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
            NW = randn(N,N);
            NB = randn(N,N);
            NZ = NW*rho.' + NB*sqrtm(eye(N)-rho*rho.');
            %NZ = rho.'*NW + sqrtm(eye(N)-rho'*rho)*NB;
            dW = NW * sqrt(dt);
            dZ = NZ * sqrt(dt);
            
            % positive semi-definitness
            if sum(eig(V_latest)>=0)<N 
                [P, Q] = eig(V_latest);
                Q = max(Q,0);
%                 V_latest = real(P*Q*P');
                V_latest = P*Q*P';
            end
            
            interet_rate(step) = max(trace(Rn * V_latest)+hn,0);
            % Update X
            sum1 = trace((An - Am) * V_latest * (An + Am));
%             sum1 = trace((An - Am) * (An + Am) * V_latest);
            mu = trace(V_latest * R) +hn-hm+  0.5 * (sum1);
            sum2 = trace((An - Am) * sqrtm(V_latest) * dW); 
            x_latest = x_latest + mu*dt + sum2;

            % Update V
            V_latest = V_latest + (beta*(sigma'*sigma) -0.5* kappa*V_latest-0.5* V_latest*kappa') * dt ...
                +0.5*(sqrtm(V_latest)* dZ*sigma +sigma'*dZ'*sqrtm(V_latest));
        
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

end
