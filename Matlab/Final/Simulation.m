%% Set parameter values for multifactor stochastic volatility model
% Recchioni, Sun, An explicitly solvable Heston model with stochastic interest rate,
% European Journal of Operational Research 249 (1), 359-377, 2016.
% Call or put parameter
% theta = -1; % 1 for call and -1 for put
Simulation_compare(1)
Simulation_compare(-1)

function Simulation_compare(theta)
% Market parameters
S0 = 1; % spot exchange rate
r_0 = [0.05,0.06]; % spot interest rates r_{i0},r_{j0}

% Contract parameters
T = 1; % maturity
%saK = 1; % strike price

% Model parameters
param_alpha = 0.5; %
d = 2; % number of volatility factors

% Volatility coefficients or weights
b_i = [0.6650 1.0985];
b_j = [1.6177 1.3588];

% Mean-reversion rate (or strength) of the volatility
kappav = [0.9418,1.7909];

% Initial volatility
v_0 = [0.1244,0.0591];

% Long-term average of the volatility
v_bar = [0.037,0.0909];

% Volatility of volatility
sigmav = [0.4912,0.08];

% Interest rate coefficients or weights
a_i = [1.004,0.000000];
a_j = [0.000000,1.006];

% Mean-reversion rate (or strength) of the interest rate
kappar = [0.02,0.02]; % lambda_i,lambda_j

% Long-term average of the interest rate
r_bar = [0.05,0.06]; % \bar{r}_i,\bar{r}_j

% Volatility of the interest rate
sigmar = [0.002,0.002]; % \eta_i,\eta_j

% Correlations
rho_v = [0.5231,-0.398];
rho_r = [-0.23,-0.81];

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

% Auxiliary parameters
b_ij_minus = b_i - b_j;
a_ij_minus = a_i - a_j;
b_ij_plus = b_i + b_j;
a_ij_plus = a_i + a_j;
b_ij_rho = rho_v.*b_ij_minus;
a_ij_rho = rho_r.*a_ij_minus;
b_ij_division = b_ij_plus./b_ij_minus;
a_ij_division = a_ij_plus./a_ij_minus;

fr = xi.^2-((alpha^2-alpha*a_ij_division')*ones(1,ngrid)-(2*alpha-a_ij_division')*1i*xi);%
fv = xi.^2-((alpha^2-alpha*b_ij_division')*ones(1,ngrid)-(2*alpha-b_ij_division')*1i*xi);%
er = (kappar.'*ones(1,ngrid)+(sigmar.*a_ij_rho)'.*(-1i*xi+alpha));
ev = (kappav.'*ones(1,ngrid)+(sigmav.*b_ij_rho)'.*(-1i*xi+alpha));
dr= (er.^2 + diag(sigmar.^2)*(diag(a_ij_minus.^2)*fr+2*diag(a_ij_division)*ones(2,1)*(alpha-1i*xi))).^0.5; %
dv= (ev.^2 + diag(sigmav.^2.*b_ij_minus.^2)*fv).^0.5;
gr = (er - dr)./ (er + dr);
gv = (ev - dv)./ (ev + dv);
CFr = exp(kappar.*r_bar./sigmar.^2*((er-dr)*T-2*log((1-gr.*exp(-dr*T))./(1-gr))) + r_0./sigmar.^2*((er-dr).*(1-exp(-dr*T))./(1-gr.*exp(-dr*T))));
CFv = exp(kappav.*v_bar./sigmav.^2*((ev-dv)*T-2*log((1-gv.*exp(-dv*T))./(1-gv))) + v_0./sigmav.^2*((ev-dv).*(1-exp(-dv*T))./(1-gv.*exp(-dv*T))));
%CF = CFr.*CFv.*exp(1i*xi*T*r_0(1));
CF = CFr.*CFv;


for K = [0.9 0.95 1 1.05 1.1]
    % Recchioni Page 6, 2016 eq. (34)
    % Sun page 31 eq. (154,155)
    factor_simple = S0*exp(-r_0(1)*T); % mixes discount and damping
    payoff = (K/S0).^(alpha+1+1i*xi)./((+1i*xi+alpha).*(+1i*xi+alpha+1));
    integrand_new = conj(payoff).*CF;
    priceS_new = factor_simple*sum(integrand_new)*dxi/(2*pi);
    if theta ==1
        fprintf('The call price of %2.2f is %4.6f\n', K,priceS_new)
    else
        fprintf('The put price of %2.2f is %4.6f\n', K,priceS_new)
    end

    % Algorithm parameters
    nblocks = 20;
    npaths = 500;
    parfor n = [4 5 6 7 8 9 10]
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
        VcMC_result = mean(VcMC);
        VpMC_result = mean(VpMC);
        scMC = sqrt(var(VcMC)/nblocks);
        spMC = sqrt(var(VpMC)/nblocks);
        
        ts = tinv([0.025  0.975],nblocks-1);      % 95 T-Score
        CI_c = VcMC_result+ ts*scMC;
        CI_p = VpMC_result+ ts*spMC;
        
        if theta ==1
            fprintf('nsteps = %3d, MCprice is %4.6f, std is %4.6f, abs err is %4.6f\n',nsteps,VcMC_result,scMC,abs(priceS_new-VcMC_result))
            if (priceS_new > CI_c(1) && priceS_new < CI_c(2))
                fprintf('Validated\n')
            else
                fprintf('Not Validated\n')
            end
        else
            fprintf('nsteps = %3d, MCprice is %4.6f, std is %4.6f, abs err is %4.6f\n',nsteps,VpMC_result,spMC,abs(priceS_new-VpMC_result))
            if (priceS_new > CI_p(1) && priceS_new < CI_p(2))
                fprintf('Validated\n')
            else
                fprintf('Not Validated\n')
            end
        end

    end
end
end