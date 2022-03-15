% Monte Carlo parameters; 
nblocks = 8000;
npaths = 5000;


%% Monte Carlo

tic;
VcMC = zeros(nblocks,1);
VpMC = zeros(nblocks,1);
for block = 1:nblocks
    VcMCb = zeros(1,npaths);
    VpMCb = zeros(1,npaths);
    for path = 1:npaths
        v_1 = zeros(nsteps,1);
        v_1 = [v_1_0;v_1]; 
        v_2 = zeros(nsteps,1);
        v_2 = [v_2_0;v_2];
        % corr (dW_v_1, dZ_v_1) = rho_v_1
        dW_v_1_1 = randn(nsteps,1);
        dW_v_1_help_1 = randn(nsteps,1);
        dZ_v_1_1 = rho_v_1 * dW_v_1_1 + (1-rho_v_1^2)^0.5 * dW_v_1_help_1;   
        % corr (dW_v_2, dZ_v_2) = rho_v_2
        dW_v_2_1 = randn(nsteps,1);
        dW_v_2_help_1 = randn(nsteps,1);
        dZ_v_2_1 = rho_v_2 * dW_v_2_1 + (1-rho_v_2^2)^0.5 * dW_v_2_help_1; 
        % Convert random numbers into Wiener processes
        dW_v_1 = dW_v_1_1 * sqrt(dt);
        dW_v_2 = dW_v_2_1 * sqrt(dt);
        dZ_v_1 = dZ_v_1_1 * sqrt(dt);
        dZ_v_2 = dZ_v_2_1 * sqrt(dt);        

        for steps = 1:nsteps
            v_1(steps+1) = max(real(v_1(steps) + chi_1 * (v_bar_1 - v_1(steps)) * dt + gamma_1 * sqrt(v_1(steps)) * dZ_v_1(steps)),0);
            v_2(steps+1) = max(real(v_2(steps) + chi_2 * (v_bar_2 - v_2(steps)) * dt + gamma_2 * sqrt(v_2(steps)) * dZ_v_2(steps)),0);
        end
        x = zeros(nsteps+1,1);

        for steps = 1:nsteps
            mu_yu = (r_i_0 - r_j_0 - 0.5 * v_1(steps) * (a_i(1)^2-a_j(1)^2) - 0.5 * v_2(steps) * (a_i(2)^2-a_j(2)^2));
             x(steps+1) = x(steps) + dt * mu_yu + v_1(steps)^0.5*(a_i(1)-a_j(1)) * dW_v_1(steps) +  v_2(steps)^0.5*(a_i(2)-a_j(2)) * dW_v_2(steps);%
        end
        S_end = S0*exp(x(end));
    
        payoffs_call = max(S_end - K,0);
        payoffs_put = max(K - S_end,0);
    
        VcMCb(1,path) = exp(-r_i_0*T)*payoffs_call;
        VpMCb(1,path) = exp(-r_i_0*T)*payoffs_put;
    end
    VcMC(block) = mean(VcMCb);
    VpMC(block) = mean(VpMCb);
end
VcMC_result = mean(VcMC);
VpMC_result = mean(VpMC);
scMC = sqrt(var(VcMC)/nblocks);
spMC = sqrt(var(VpMC)/nblocks);

cputime_MC = toc;
fprintf('%22s%14.10f%14.10f%14.3f\n','Monte Carlo',VcMC_result,VpMC_result,cputime_MC)
fprintf('%22s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)