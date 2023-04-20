function [simulated_call, simulated_put] = europeanSimulation(params,european,n,S0,T,K)
% Model parameters
a_i = diag(params.a_i);
a_j = diag(params.a_j);
kappa = diag(params.kappa); 
y_bar = diag(params.y_bar); 
sigma = diag(params.sigma); 
y_0 = diag(params.y_0);
rho = diag(params.rho);
hm = diag(params.hm);
hn = diag(params.hn);
h = hm-hn;
d = 4;
% Number of simulations 
nblocks = 50;
npaths = 100;
% Number of steps 
nsteps = 50*n;
dt = T/nsteps;
% Initialize calls
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    for path = 1:npaths
        yPaths = y_0(:,:,ones(1,nsteps+1));
        % corr (dW_v, dZ_v) = rho_v
        dW_1 = randn(nsteps,d);
        dW_2 = randn(nsteps,d);
        dW_3 = dW_1*rho + dW_2*sqrt(eye(d)-rho.^2);

        dW = dW_1 * sqrt(dt);
        dZ = dW_3 * sqrt(dt);

        % dv = chi * (v_bar - v) * dt + gamma * \sqrt(v) * dZ_v
        for step = 1:nsteps
            yPaths(:,:,step+1) = max(yPaths(:,:,step) + (y_bar - yPaths(:,:,step)) * kappa * dt ...
                + yPaths(:,:,step).^0.5 * sigma * diag(dZ(step,:)),0);
        end

        sum_1 = zeros(nsteps,1);
        sum_2 = zeros(nsteps,1);
        mu = zeros(nsteps,1);
        x = zeros(nsteps+1,1);
        % Valuation for dx
        for step = 1:nsteps
            sum_1(step) = trace(yPaths(:,:,step) * (a_i + a_j)' * (a_i - a_j) );
            mu(step) = trace(yPaths(:,:,step) * h) +  0.5 * (sum_1(step));
            sum_2(step) = trace(yPaths(:,:,step).^0.5 * (a_i-a_j)' * diag(dW(step,:)) );
            x(step+1) = x(step) + mu(step)*dt + sum_2(step);
        end

        % Calculate the Smax and Smin

        S_end = S0*exp(x(end));

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);

        % Countinuous Compounding 
        r_sum = 0;
        for step = 1:nsteps
            r_sum = r_sum + trace(yPaths(:,:,step) * hm) *dt;
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

ts = tinv([0.025  0.975],nblocks-1);      % 95 T-Score
CI_c = simulated_call+ ts*scMC;
CI_p = simulated_put+ ts*spMC;


fprintf('nsteps = %3d, MCprice is %4.6f, std is %4.6f, abs err is %4.6f  ',nsteps,simulated_call,scMC,abs(european.call_price-simulated_call))
if (european.call_price > CI_c(1) && european.call_price < CI_c(2))
    fprintf('Validated\n')
else
    fprintf('Not Validated\n')
end

fprintf('nsteps = %3d, MCprice is %4.6f, std is %4.6f, abs err is %4.6f  ',nsteps,simulated_put,spMC,abs(european.put_price-simulated_put))
if (european.put_price > CI_p(1) && european.put_price < CI_p(2))
    fprintf('Validated\n')
else
    fprintf('Not Validated\n')
end

end