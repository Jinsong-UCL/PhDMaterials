function [simulated_call, simulated_put] = europeanSimulation(params,european,n,S0,T,K,r_0)

a_i = params.a_i;
a_j = params.a_j;
kappar = params.kappar;
r_bar = params.r_bar;
sigmar = params.sigmar;
b_i = params.b_i;
b_j = params.b_j;
kappav = params.kappav;
v_0 = params.v_0;
v_bar = params.v_bar;
sigmav = params.sigmav;
rho_v = params.rho_v;
rho_r = params.rho_r;
d = 2;
param_alpha = 0.5;

nblocks = 20;
npaths = 50;

nsteps = 2^n;
dt = T/nsteps;
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
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
            v(steps+1,:) = max(v(steps,:) + (v_bar - v(steps,:)) * diag(kappav) * dt ...
                + sqrt(v(steps,:)) .* dZ_v(steps,:)* diag(sigmav),0);
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
            r(steps+1,:) = max(r(steps,:) + (r_bar - r(steps,:)) * diag(kappar) * dt ...
                + r(steps,:).^param_alpha .* dZ_r(steps,:)* diag(sigmar),0);
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
            sum_v_1(steps) = v(steps,:) * (b_i.^2' - b_j.^2');
            % sum_r_1(steps) = b_i^2 * r_i(steps)^(2*param_alpha) - b_j^2 * r_j(steps)^(2*param_alpha);
            sum_r_1(steps) = r(steps,:) * (a_i.^2' - a_j.^2');

            mu(steps) = r(steps,1) - r(steps,2) + 0.5 * (sum_v_1(steps) + sum_r_1(steps));

            % Block for v
            sum_v_2(steps) = v(steps,:).^0.5 .* dW_v(steps,:) * (b_i'-b_j');

            % Block for r
            sum_r_2(steps) = r(steps,:).^param_alpha .* dW_r(steps,:)*(a_i'-a_j');

            % dx = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i *
            % diag_r^alpha * dW_r ???????????????
            x(steps+1) = x(steps) + mu(steps)*dt + sum_v_2(steps) + sum_r_2(steps);
        end

        % Calculate the Smax and Smin

        S_end = S0*exp(x(end));

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);

        VcMCb(1,path) = exp(-r_0(1)*T)*payoffs_call;
        VpMCb(1,path) = exp(-r_0(1)*T)*payoffs_put;

        % MC_result(1:nsteps+1,path) = S0*exp(x_i_j);
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