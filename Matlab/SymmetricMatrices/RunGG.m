%% Warning option
[msgStr,warnId] = lastwarn;
warnStruct = warning('off',warnId);

%% Market parameters 
marketStruct.S0 = 1;
marketStruct.K = 1;
marketStruct.d = 2;
marketStruct.T = 1;
% The following numbers are from Gnoatto and Grasselli 2014
%% Diagonal parameters
% paramStruct.beta = 3.1442;
% paramStruct.am = [0.7764 0.0;0.0 0.9639];
% paramStruct.an = [0.6679 0.0;0.0 0.8520];
% paramStruct.hm = [0.2725 0.0;0.0 0.4726];
% paramStruct.hn = [0.1841 0.0;0.0 0.4700];
% paramStruct.y_0 = [0.1688 0.0;0.0 0.3169];
% paramStruct.rho = [-0.5417 0.0;0.0 -0.4834];
% paramStruct.kappa = [1.0426,0.0;0.0,0.8778]; 
% paramStruct.sigma = [0.4364,0.0;0.0,0.7362];
%% Symmetric parameters
paramStruct.beta = 3.1442;
paramStruct.am = [0.7764 0.4837;0.4837 0.9639];
paramStruct.an = [0.6679 0.6277;0.6277 0.8520];
paramStruct.hm = [0.2725 0.0804;0.0804 0.4726];
paramStruct.hn = [0.1841 0.0155;0.0155 0.4761];
paramStruct.y_0 = [0.1688 0.1708;0.1708 0.3169];
paramStruct.rho = [-0.5417 0.1899;0.1899 -0.4834];
paramStruct.kappa = [1.0426,0.6764;0.6764,0.8778]; 
paramStruct.sigma = [0.4364,0.1914;0.1914,0.7362];

% rho test this is to test if 
if sum(eig(eye(marketStruct.d)-paramStruct.rho*paramStruct.rho.')>=0) == marketStruct.d 
    fprintf("rho is valid, please continue.\n");
else
    fprintf("rho is not valid, please check rho first.\n");
end

%% Simulation without constant dispalcements in the interets rates
%[simulated_call, simulated_put, phi_empirical] = GGsimulation(marketStruct,paramStruct);

%% CF
fprintf("This is the result of CF in our format\n")
[call_price_H] = HCF(marketStruct,paramStruct,1);
[put_price_H] = HCF(marketStruct,paramStruct,-1);

%% CFG
fprintf("This is the result of CF in Da Fonseca's format\n")
[call_price] = GCF(marketStruct,paramStruct,1);
[put_price] = GCF(marketStruct,paramStruct,-1);

%% ECF and ACF
%[statement] = ECFACF(marketStruct,paramStruct,phi_empirical);

