% Market parameters
S0 = 1; % spot exchange rate
r_0 = [0.05,0.06]; % spot interest rates r_{i0},r_{j0}

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Interest rate coefficients or weights
parStruct.a_i = [1.004,0.000000];
parStruct.a_j = [0.000000,1.006];

% Mean-reversion rate (or strength) of the interest rate
parStruct.kappar = [0.02,0.02]; % lambda_i,lambda_j

% Long-term average of the interest rate
parStruct.r_bar = [0.05,0.06]; % \bar{r}_i,\bar{r}_j

% Volatility of the interest rate
parStruct.sigmar = [0.002,0.002]; % \eta_i,\eta_j

% Volatility coefficients or weights
parStruct.b_i = [0.6650 1.0985];
parStruct.b_j = [1.6177 1.3588];

% Mean-reversion rate (or strength) of the volatility
parStruct.kappav = [0.9418,1.7909];
% Initial volatility
parStruct.v_0 = [0.1244,0.0591];
% Long-term average of the volatility
parStruct.v_bar = [0.037,0.0909];

% Volatility of volatility
parStruct.sigmav = [0.4912,0.08];

% Correlations
parStruct.rho_v = [0.5231,-0.398];
parStruct.rho_r = [-0.23,-0.81];

%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% European Options
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
european.call_price = europeanPricing(parStruct,1,S0,T,K,r_0);
european.put_price = europeanPricing(parStruct,-1,S0,T,K,r_0);


%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Barrier Options
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
barrier.call_price = barrierPricing(parStruct,exp(-9),exp(9),1,S0,T,K,r_0);
barrier.put_price = barrierPricing(parStruct,exp(-9),exp(9),-1,S0,T,K,r_0);

%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Simulations
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Number of steps 2^n
parfor n = [4 5 6 7 8 9 10]
    [simulated_call, simulated_put] = europeanSimulation(parStruct,european,n,S0,T,K,r_0);
end

