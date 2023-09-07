function [simulated_call, simulated_put] = europeanSimulation(params,european,n,S0,T,K)
% Model parameters
a_i = params.a_i;
a_j = params.a_j;
kappa = params.kappa; 
y_bar = params.y_bar; 
sigma = params.sigma; 
y_0 = params.y_0;
rho = params.rho;
hm = params.hm;
hn = params.hn;
h = hm-hn;
d = 4;
% Number of simulations 
nblocks = 100;
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
        % Increments of the arithmetic Brownian motion X(t) = log(S(t)/S(0))
        % dX = muABM*dt + sigma*sqrt(dt)*randn(ndates,nsample);
        % dS/S = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r

        % Generate a random path using the above model
        % Volatility part
        % Model variables
        y = zeros(nsteps+1,d);
        y(1,:) = y_0;
        % corr (dW_v, dZ_v) = rho_v
        dW_1 = randn(nsteps,d);
        dW_2 = randn(nsteps,d);
        dW_3 = dW_1*diag(rho) + dW_2*diag((1-rho.^2).^0.5);

        dW = dW_1 * sqrt(dt);
        dZ = dW_3 * sqrt(dt);

        % dv = chi * (v_bar - v) * dt + gamma * \sqrt(v) * dZ_v
        for step = 1:nsteps
            y(step+1,:) = max(y(step,:) + (y_bar - y(step,:)) * diag(kappa) * dt ...
                + sqrt(y(step,:)) .* dZ(step,:)* diag(sigma),0);
        end

        sum_y_1 = zeros(nsteps,1);
        sum_y_2 = zeros(nsteps,1);
        mu = zeros(nsteps,1);
        x = zeros(nsteps+1,1);
        % Valuation for dx
        for step = 1:nsteps
            sum_y_1(step) = y(step,:) * (a_i.^2' - a_j.^2'); % plus
            %sum_y_1(step) = y(step,:) * (a_i' - a_j').^2; % minus
            mu(step) = sum(y(step,:).*h) + 0.5 * (sum_y_1(step));
            sum_y_2(step) = y(step,:).^0.5 .* dW(step,:) * (a_i'-a_j');
            x(step+1) = x(step) + mu(step)*dt + sum_y_2(step);
        end

        % Calculate the Smax and Smin

        S_end = S0*exp(x(end));

        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);

        % Countinuous Compounding 
        r_sum = 0;
        for step = 1:nsteps
            r_sum = r_sum + trace(y(step,:) * hm') *dt;
        end
        dfactor = exp(- r_sum);
        %fprintf("The difference is %f", abs(dfactor2-dfactor))

        VcMCb(1,path) = dfactor*payoffs_call;
        VpMCb(1,path) = dfactor*payoffs_put;
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
% fprintf('&& %3d & %4.6f & %4.6f & %4.6f  ',nsteps,simulated_call,scMC,abs(european.call_price-simulated_call))
if (european.call_price > CI_c(1) && european.call_price < CI_c(2))
    %fprintf('Validated\n')
    fprintf('\\\\ \n')
else
    %fprintf('Not Validated\n')
    fprintf('*\\\\ \n')
end

fprintf('nsteps = %3d, MCprice is %4.6f, std is %4.6f, abs err is %4.6f  ',nsteps,simulated_put,spMC,abs(european.put_price-simulated_put))
% fprintf('&& %3d & %4.6f & %4.6f & %4.6f  ',nsteps,simulated_put,spMC,abs(european.put_price-simulated_put))
if (european.put_price > CI_p(1) && european.put_price < CI_p(2))
    %fprintf('Validated\n')
    fprintf('\\\\ \n')
else
    %fprintf('Not Validated\n')
    fprintf('*\\\\ \n')
end
% figure(1)
% plot(xi,real(CF_E),xi,real(GCF_A),"--")
% axis([-20 20 -0.1 0.7])
% title('Real part of the discounted characteristic function for the matrix Heston model')
% xlabel('\xi')
% legend('Empirical','Analytical')
% savefig("Real part of the discounted characteristic function for the matrix Heston model(asymmetric).fig")
% 
% figure(2)
% plot(xi,imag(CF_E),xi,imag(GCF_A),"--")
% axis([-20 20 -0.3 0.3])
% title('Imaginary part of the discounted characteristic function for the matrix Heston model')
% xlabel('\xi')
% legend('Empirical','Analytical')
% savefig("Imaginary part of the discounted characteristic function for the matrix Heston model(asymmetric).fig")

end