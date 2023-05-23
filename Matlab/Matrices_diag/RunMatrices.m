%% Warning option
[msgStr,warnId] = lastwarn;
warnStruct = warning('off',warnId);

%% Market parameters 
marketStruct.S0 = 1;
K = 1;
marketStruct.d = 2;
marketStruct.T = 1;
% The following numbers are from Gnoatto and Grasselli 2014
%% Symmetric parameters
paramStruct.beta = 3.144;
paramStruct.An = [0.7764 0;0 0.9639];
paramStruct.Am = [0.6679 0;0 0.8520];
paramStruct.Rn = [0.2725 0;0 0.4726];
paramStruct.Rm = [0.1841 0;0 0.4761];
paramStruct.V_0 = [0.1688 0;0 0.3169];
paramStruct.rho = [-0.5417 0;0 -0.4834];
paramStruct.kappa = [1.0426,0;0,0.8778]; 
paramStruct.sigma = [0.4368,0;0,0.7362];
%paramStruct.sigma = [0.4368,0.1914;0.1914,0.7362];
%paramStruct.sigma = [0.4368,0.4214;0.1914,0.7362];
%paramStruct.sigma = [0.441736679745699,0.355320710364647;0.355320710364647,0.749833853077106];
% Number of simulations
paramStruct.nblocks = 100;
paramStruct.npaths =  100;
%% Fourier parameters
fourierStruct.xwidth = 20; % width of the support in real space
fourierStruct.ngrid = 2^12; % number of grid points

% Grids in real and Fourier space
fourierStruct.B = fourierStruct.xwidth/2; % upper bound of the support in real space
fourierStruct.dxi = pi/fourierStruct.B; % Nyquist relation: grid step in Fourier space
fourierStruct.xi = fourierStruct.dxi*(-fourierStruct.ngrid/2:fourierStruct.ngrid/2-1); % grid in Fourier space

% rho test this is to test if 
if sum(eig(eye(marketStruct.d)-paramStruct.rho*paramStruct.rho.')>=0) == marketStruct.d 
    fprintf("rho is valid.\n");
else
    error("rho is not valid, please check rho first.\n");
end
 
%% Simulation without constant dispalcements in the interets rates
[simulated_call, simulated_put,~,~,CF_E] = GGsimulation(marketStruct,paramStruct,fourierStruct,K,3);

%% CF
fprintf("This is the result of GHCF\n")
call_price_H = HCF(marketStruct,paramStruct,fourierStruct,K,1);
put_price_H = HCF(marketStruct,paramStruct,fourierStruct,K,-1);

%% CFG
fprintf("This is the result of GGCF\n")
call_price_G = GCF(marketStruct,paramStruct,fourierStruct,K,1);
put_price_G = GCF(marketStruct,paramStruct,fourierStruct,K,-1);

%% ECF and ACF
[statement] = ECFACF(marketStruct,paramStruct,fourierStruct,CF_E);

