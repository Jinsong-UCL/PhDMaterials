import torch

# Market parameters
S0 = 1; # spot exchange rate
r_0 = torch.FloatTensor([0.02,0.01]); # spot interest rates r_{i0},r_{j0}

# Contract parameters
T = 1; # maturity
K = 1; # strike price

# Model parameters

param_alpha = 0.5;
d = 2; # number of volatility factors

# Volatility coefficients or weights
a_i =torch.FloatTensor( [0.6650, 1.0985]);
a_j = torch.FloatTensor([1.6177, 1.3588]);

# Mean-reversion rate (or strength) of the volatility
chi = torch.FloatTensor([0.9418,1.7909]);

# Initial volatility
v_0 = torch.FloatTensor([0.1244,0.0391]);

# Long-term average of the volatility
v_bar = torch.FloatTensor([0.037,0.0909]);

# Volatility of volatility
gamma = torch.FloatTensor([0.4912,1]);

# Interest rate coefficients or weights
b_i = torch.FloatTensor([1.0000004,0.000000]);
b_j = torch.FloatTensor([0.000000,0.0000006]);

# Mean-reversion rate (or strength) of the interest rate
Lambda = torch.FloatTensor([0.01,0.02]); # lambda_i,lambda_j

# Long-term average of the interest rate
r_bar = torch.FloatTensor([0.02,0.01]); # \bar{r}_i,\bar{r}_j

# Volatility of the interest rate
eta = torch.FloatTensor([0.001,0.002]); # \eta_i,\eta_j

# Correlations
rho_v = torch.FloatTensor([-0.5231,-0.398]);
rho_r = torch.FloatTensor([-0.23,-0.81]);


# Algorithm parameters
nsteps = 100;
nblocks = 200;
npaths = 1000;


# Monte Carlo
dt = T/nsteps;
VcMC = torch.zeros(nblocks,1);
VpMC = torch.zeros(nblocks,1);
for block in range (nblocks):
    MC_result = torch.zeros(nsteps+1,npaths);
    VcMCb = torch.zeros(1,npaths);
    VpMCb = torch.zeros(1,npaths);
    for path in range(npaths):
        # Increments of the arithmetic Brownian motion X(t) = log(S(t)/S(0))
        # dX = muABM*dt + sigma*sqrt(dt)*randn(ndates,nsample);
        # dS/S = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r

        # Generate a random path using the above model
        # Volatility part
        # Model variables
        v = torch.zeros(nsteps+1,d);
        v[1,] = v_0;
        # corr (dW_v, dZ_v) = rho_v
        dW_v_1 = torch.randn(nsteps,d);
        dW_v_2 = torch.randn(nsteps,d);
        dW_v_3 = dW_v_1*torch.diag(rho_v) + dW_v_2*torch.diag((1-rho_v^2)^0.5);

        dW_v = dW_v_1 * torch.sqrt(dt);
        dZ_v = dW_v_3 * torch.sqrt(dt);

        # dv = chi * (v_bar - v) * dt + gamma * \sqrt(v) * dZ_v
        for steps in range(nsteps):
            v[steps+1,] = max(v[steps,] + (v_bar - v[steps,]) * torch.diag(chi) * dt
                              + torch.sqrt(v[steps,]) * dZ_v[steps,]* torch.diag(gamma),0);

        # Interest rate part
        # Model variables
        r = torch.zeros(nsteps+1,2);
        r[1,] = r_0;

        # corr (dW_r_i, dZ_r_i) = rho_r_i
        dW_r_1 = torch.randn(nsteps,2);
        dW_r_2 = torch.randn(nsteps,2);
        dW_r_3 = dW_r_1*torch.diag(rho_r) + dW_r_2*torch.diag((1-rho_r^2)^0.5);

        dW_r = dW_r_1 * torch.sqrt(dt);
        dZ_r = dW_r_3 * torch.sqrt(dt);

        # dr = lambda * (r_bar - r) * dt + eta * r^alpha * dZ_r
        for steps in range(nsteps):
            r[steps+1,] = max(r[steps,] + (r_bar - r[steps,]) * torch.diag(Lambda) * dt
                              + r[steps,]^param_alpha * dZ_r[steps,]* torch.diag(eta),0);

        sum_v_1 = torch.zeros(nsteps,1);
        sum_v_2 = torch.zeros(nsteps,1);
        sum_r_1 = torch.zeros(nsteps,1);
        sum_r_2 = torch.zeros(nsteps,1);
        mu = torch.zeros(nsteps,1);
        x = torch.zeros(nsteps+1,1);

        # Valuation for dx
        for steps in range(nsteps):
            # Block for MuYu
            # sum_v_1(steps) = (a_i(1)^2 - a_j(1)^2) * v_1(steps) + (a_i(2)^2 - a_j(2)^2) * v_2(steps);
            sum_v_1[steps] = v[steps,] * (a_i^2 - a_j^2);
            # sum_r_1(steps) = b_i^2 * r_i(steps)^(2*param_alpha) - b_j^2 * r_j(steps)^(2*param_alpha);
            sum_r_1[steps] = r[steps,] * (b_i^2 - b_j^2);

            mu[steps] = r[steps,1] - r[steps,2] + 0.5 * (sum_v_1[steps] + sum_r_1[steps]);

            # Block for v
            sum_v_2[steps] = torch.dot(v[steps,:]^0.5 * dW_v[steps,] , (a_i-a_j));###

            # Block for r
            sum_r_2[steps] = torch.dot(r[steps,]^param_alpha * dW_r[steps,],(b_i-b_j));###

            # dx = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r
            x[steps+1] = x[steps] + mu[steps]*dt + sum_v_2[steps] + sum_r_2[steps];

        # Calculate the Smax and Smin

        S_end = S0*torch.exp(x(-1));

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);

        VcMCb[1,path] = torch.exp(-r_0(1)*T)*payoffs_call;
        VpMCb[1,path] = torch.exp(-r_0(1)*T)*payoffs_put;

       # MC_result(1:nsteps+1,path) = S0*exp(x_i_j);

    VcMC[block] = torch.mean(VcMCb);
    VpMC[block] = torch.mean(VpMCb);
    #fprintf('%14.10f\n',block);

VcMC_result = torch.mean(VcMC);
VpMC_result = torch.mean(VpMC);
scMC = torch.sqrt(torch.var(VcMC)/nblocks);
spMC = torch.sqrt(torch.var(VpMC)/nblocks);