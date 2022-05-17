function [payoffs_call,payoffs_put] = MC_GPU(nsteps,d,v_0,dt,rho_v,r_0,rho_r,a_i,a_j,b_i,b_j,param_alpha,S0,K)
        % Increments of the arithmetic Brownian motion X(t) = log(S(t)/S(0))
        % dX = muABM*dt + sigma*sqrt(dt)*randn(ndates,nsample);
        % dS/S = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r    
    
        % Generate a random path using the above model
        % Volatility part
        % Model variables
        v = zeros(nsteps+1,d);
        v(1,:) = v_0;
        % corr (dW_v, dZ_v) = rho_v
        dW_v_1 = randn(nsteps,d,"gpuArray");
        dW_v_2 = randn(nsteps,d,"gpuArray");
        dW_v_3 = dW_v_1*diag(rho_v) + dW_v_2*diag((1-rho_v.^2).^0.5);

        dW_v = dW_v_1 * sqrt(dt);
        dZ_v = dW_v_3 * sqrt(dt);
        
        % dv = chi * (v_bar - v) * dt + gamma * \sqrt(v) * dZ_v
        for steps = 1:nsteps
            v(steps+1,:) = max(v(steps,:) + (v_bar - v(steps,:)) * diag(chi) * dt ...
                + sqrt(v(steps,:)) .* dZ_v(steps,:)* diag(gamma),0);
        end
    
        % Interest rate part
        % Model variables
        r = zeros(nsteps+1,2);
        r(1,:) = r_0;

        % corr (dW_r_i, dZ_r_i) = rho_r_i
        dW_r_1 = randn(nsteps,2,"gpuArray");
        dW_r_2 = randn(nsteps,2,"gpuArray");
        dW_r_3 = dW_r_1*diag(rho_r) + dW_r_2*diag((1-rho_r.^2).^0.5);

        dW_r = dW_r_1 * dt;
        dZ_r = dW_r_3 * dt;

        % dr = lambda * (r_bar - r) * dt + eta * r^alpha * dZ_r
        for steps = 1:nsteps   
            r(steps+1,:) = max(r(steps,:) + (r_bar - r(steps,:)) * diag(lambda) * dt ...
                + r(steps,:).^param_alpha .* dZ_r(steps,:)* diag(eta),0);
        end
        
        sum_v_1 = zeros("like",nsteps,1);
        sum_v_2 = zeros("like",nsteps,1);
        sum_r_1 = zeros("like",nsteps,1);
        sum_r_2 = zeros("like",nsteps,1);
        mu = zeros("like",nsteps,1);
        x = zeros("like",nsteps+1,1);
        % Valuation for dx
        for steps = 1:nsteps
            % Block for MuYu
            % sum_v_1(steps) = (a_i(1)^2 - a_j(1)^2) * v_1(steps) + (a_i(2)^2 - a_j(2)^2) * v_2(steps);
            sum_v_1(steps) = v(steps,:) * (a_i.^2' - a_j.^2');
            % sum_r_1(steps) = b_i^2 * r_i(steps)^(2*param_alpha) - b_j^2 * r_j(steps)^(2*param_alpha);
            sum_r_1(steps) = r(steps,:) * (b_i.^2' - b_j.^2');

            mu(steps) = r(steps,1) - r(steps,2) + 0.5 * (sum_v_1(steps) + sum_r_1(steps));
            % muYu(steps,1) = r_i(steps,1) - r_j(steps,1) - 0.5 * (sum_v_1(steps) + sum_r_1(steps));
            
            % Block for v
            sum_v_2(steps) = v(steps,:).^0.5 .* dW_v(steps,:) * (a_i'-a_j');
            
            % Block for r
            sum_r_2(steps) = r(steps,:).^param_alpha .* dW_r(steps,:)*(b_i'-b_j');
    
            % dx = (r_0 - r_i)*dt - a_i * diag_v^0.5 * dW_v - b_i * diag_r^alpha * dW_r 
            x(steps+1) = x(steps) + mu(steps)*dt + sum_v_2(steps) + sum_r_2(steps);
        end
    
        % Calculate the Smax and Smin
    
        S_end = S0*exp(x(end));
    
        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);
end

