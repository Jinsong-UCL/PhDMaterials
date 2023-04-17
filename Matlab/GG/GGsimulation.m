% Simulation based on our notations 
S0 = 1;
K = 1;
% The following numbers are from Gnoatto and Grasselli 2014
%am = [0.7764 0.4837;0.4837 0.9639];
%an = [0.6679 0.6277;0.6277 0.8520];
%hm = [0.2725 0.0804;0.0804 0.4726];
%hn = [0.1841 0.0155;0.0155 0.4761];
%y_0 = [0.1688 0.1708;0.1708 0.3169];
am = [0.7764 0.0;0.0 0.9639];
an = [0.6679 0.0;0.0 0.8520];
hm = [0.2725 0.0;0.0 0.4726];
hn = [0.1841 0.0;0.0 0.4700];
y_0 = [0.1688 0.0;0.0 0.3169];
%y_bar =
%rho = [-0.5417 0.1899;-0.1170 -0.4834];
%kappa = [1.0426,0.6764;0.9880,0.8778]; 
%sigma = [0.4364,0.1914;0.4966,0.7362];

%rho = [-0.5417 0.1899;0.1899 -0.4834];
%kappa = [1.0426,0.6764;0.6764,0.8778]; 
%sigma = [0.4364,0.1914;0.1914,0.7362];

rho = [-0.5417 0.0;0.0 -0.4834];
kappa = [1.0426,0.0;0.0,0.8778]; 
sigma = [0.4364,0.0;0.0,0.7362];

h = hm - hn;
beta = 3.1442;
T = 1;
d = 2;
% Number of simulations 
nblocks = 100;
npaths = 10;
% Number of steps 
nsteps = 100;
dt = T/nsteps;
% Initialize calls
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
tic;
[msgStr,warnId] = lastwarn;
warnStruct = warning('on',warnId);

for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    for path = 1:npaths
        yPaths = y_0(:,:,ones(1,nsteps+1));
        % W and Z are matrix Brownian motions
        dW_1 = randn(d,d,nsteps);
        dW_2 = randn(d,d,nsteps);
        dW_3 = zeros(d,d,nsteps);
        for step = 1:nsteps
            dW_3(:,:,step) = dW_1(:,:,step)*rho + dW_2(:,:,step)*sqrtm(eye(d)-rho*rho.');
        end

        dW = dW_1 * sqrt(dt);
        dZ = dW_3 * sqrt(dt);

        % dv = chi * (v_bar - v) * dt + gamma * \sqrt(v) * dZ_v
        for step = 1:nsteps
            yPaths(:,:,step+1) = max(yPaths(:,:,step) + (0.25*beta*(sigma*sigma') + kappa*yPaths(:,:,step) + yPaths(:,:,step)*kappa) * dt ...
                +0.5*(sigma*sqrtm(yPaths(:,:,step))* dZ(:,:,step) + dZ(:,:,step)'*sqrtm(yPaths(:,:,step))'*sigma'),0);
        end

        sum_1 = zeros(nsteps,1);
        sum_2 = zeros(nsteps,1);
        mu = zeros(nsteps,1);
        x = zeros(nsteps+1,1);
        % Valuation for dx
        for step = 1:nsteps
            sum_1(step) = trace((am - an) * yPaths(:,:,step)*  (am + an)');
            mu(step) = trace(yPaths(:,:,step) * h) +  0.5 * (sum_1(step));
            sum_2(step) = trace((am - an) *sqrtm(yPaths(:,:,step)) *  dW(:,:,step) );
            x(step+1) = x(step) + mu(step)*dt + sum_2(step);
        end

        % Calculate the Smax and Smin

        S_end = S0*exp(x(end));

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);

        % Countinuous Compounding 
        r_sum = 0;
        for step = 1:nsteps
            r_sum = r_sum + (trace(yPaths(:,:,step) * hm)) *dt;
        end
        dfactor = exp(- r_sum);
        
        VcMCb(1,path) = dfactor*payoffs_call;
        VpMCb(1,path) = dfactor*payoffs_put;
    end
    VcMC(block) = mean(VcMCb);
    VpMC(block) = mean(VpMCb);
end
simulated_call = mean(VcMC);
simulated_put = mean(VpMC);
scMC = sqrt(var(VcMC)/nblocks);
spMC = sqrt(var(VpMC)/nblocks);

cputime_MC = toc;
fprintf('%20s%14.10f%14.10f%14.10f\n','Monte Carlo',simulated_call,simulated_put,cputime_MC)
fprintf('%20s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)