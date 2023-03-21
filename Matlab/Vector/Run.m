% Market parameters
S0 = 1; % spot exchange rate
r_0 = [0.05,0.06]; % spot interest rates r_{i0},r_{j0}

% Contract parameters
T = 1; % maturity
K = 1; % strike price

% Interest rate coefficients or weights
a_i1 = [1.004,000000];
a_j1 = [0000000,1.0006];


% Mean-reversion rate (or strength) of the interest rate
kappar = [0.02,0.02]; % lambda_i,lambda_j

% Long-term average of the interest rate
r_bar = [0.05,0.06]; % \bar{r}_i,\bar{r}_j

% Volatility of the interest rate
sigmar = [0.002,0.002]; % \eta_i,\eta_j

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

% Correlations
rho_v = [0.5231,-0.398];
rho_r = [-0.23,-0.81];


parStruct.a_i = [a_i1,b_i];
parStruct.a_j = [a_j1,b_j];
parStruct.kappa = [kappar,kappav];
parStruct.sigma = [sigmar,sigmav];
parStruct.y_0 = [r_0, v_0];
parStruct.y_bar = [r_bar, v_bar];
parStruct.rho = [rho_r,rho_v];
parStruct.hm = [1,0,0,0];
parStruct.hn = [0,1,0,0];

%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% European Options
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fprintf("The Euripean option prices are: \n")
european.call_price = europeanPricing(parStruct,1,S0,T,K,r_0);
european.put_price = europeanPricing(parStruct,-1,S0,T,K,r_0);


%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Barrier Options
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fprintf("The Barrier option prices are: \n")
barrier.call_price = barrierPricing(parStruct,exp(-9),exp(9),1,S0,T,K,r_0);
barrier.put_price = barrierPricing(parStruct,exp(-9),exp(9),-1,S0,T,K,r_0);



%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Simulations
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[simulated_call, simulated_put] = europeanSimulation(parStruct,european,1,S0,T,K,r_0);

% Number of steps 2^n
%parfor n = [4 5 6 7]
%    [simulated_call, simulated_put] = europeanSimulation(parStruct,european,n,S0,T,K,r_0);
%end

