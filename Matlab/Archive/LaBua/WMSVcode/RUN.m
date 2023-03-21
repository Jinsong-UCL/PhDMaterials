%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Calibration and advanced simulation schemes for the Wishart Stochastic Volatility model
% (c) G. La Bua, D. Marazzina
% Corresponding Author: daniele.marazzina@polimi.it
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clear

%% Step 1: define parameters
% Model parameters
parStruct.Sigma0 = [0.0794 0.0038;0.0038 0.0003];
parStruct.M = [-0.7020 0.0893;0.0893 -0.9895];
parStruct.Q = [0.2703 -0.0198;0.0317 0.0879];
parStruct.R = [-0.7056 -0.0090;-0.0277 -0.5293];
parStruct.beta = 1.0405;
parStruct.r = 0; parStruct.q = 0;
t = 1; S0 = 100; K=S0;
% Discretization parameters
nSims = 1e4; StepsVec = 50;
z95 = norminv(1-0.05/2);

%% Step 2: pricing
[v0,parH(1),parH(2),parH(3),rho] = paramsWMSVtoHeston(t,parStruct.Sigma0,parStruct.M,parStruct.Q,parStruct.R,parStruct.beta,2);
SPaths_SIA = Paths_WMSV_SIA_MC(parStruct,parH, S0, t, StepsVec, nSims);
SPaths_GVA = Paths_WMSV_GVA_MC(parStruct, S0, t, StepsVec, nSims);

Price_SIA = exp(-parStruct.r*StepsVec)*max(SPaths_SIA(end,:)-K,0);
stdError_SIA= sqrt(var(Price_SIA)/nSims);
Price_SIA = mean(Price_SIA)
CI_SIA = [Price_SIA - z95*stdError_SIA, Price_SIA + z95*stdError_SIA]

Price_GVA = exp(-parStruct.r*StepsVec)*max(SPaths_GVA(end,:)-K,0);
stdError_GVA = sqrt(var(Price_GVA)/nSims);
Price_GVA = mean(Price_GVA)
CI_GVA = [Price_GVA - z95*stdError_GVA, Price_GVA + z95*stdError_GVA]
       
