%% Warning option
[msgStr,warnId] = lastwarn;
warnStruct = warning('off',warnId);

%% Ramdon seed option
rng(0530)

%% Write output to files
Fid = fopen("result\sigmaTest.txt","a");
%% Market parameters 
marketStruct.S0 = 1;
marketStruct.d = 2;
marketStruct.T = 1;
% The following numbers are from Gnoatto and Grasselli 2014
%% Model parameters
paramStruct.beta = 3.144;
paramStruct.An = [0.7764 0.4837;0.4837 0.9639];
paramStruct.Am = [0.6679 0.6277;0.6277 0.8520];
paramStruct.Rn = [0.2725 0.0804;0.0804 0.4726];
paramStruct.Rm = [0.1841 0.0155;0.0155 0.4761];
paramStruct.V_0 = [0.1688 0.1708;0.1708 0.3169];
paramStruct.rho = [-0.5417 0.1899;-0.1170 -0.4834];
paramStruct.kappa = [1.0426,0.6764;0.9880,0.8778]; 
paramStruct.sigma = [0.4368,0.2014;0.3514,0.7362];
% Number of simulations
paramStruct.nblocks = 100;
paramStruct.npaths =  300;

%% Fourier parameters
fourierStruct.xwidth = 20; % width of the support in real space
fourierStruct.ngrid = 2^12; % number of grid points

% Grids in real and Fourier space
fourierStruct.N = fourierStruct.ngrid/2;
fourierStruct.B = fourierStruct.xwidth/2; % upper bound of the support in real space
fourierStruct.dxi = pi/fourierStruct.B; % Nyquist relation: grid step in Fourier space
fourierStruct.xi = fourierStruct.dxi*(-fourierStruct.N:fourierStruct.N-1); % grid in Fourier space

% rho test this is to test if 
if sum(eig(eye(marketStruct.d)-paramStruct.rho*paramStruct.rho.')>=0) == marketStruct.d 
    fprintf("rho is valid.\n");
else
    error("rho is not valid, please check rho first.\n");
end


for K = [0.9 0.95 1 1.05 1.1]
    call_price = HCF(marketStruct,paramStruct,fourierStruct,K,1);
    put_price = HCF(marketStruct,paramStruct,fourierStruct,K,-1);

    parfor n = [4 5 6 7]
        warnStruct = warning('off','all');
        [simulated_call, simulated_put,scMC,spMC, ~] = GGsimulation(marketStruct,paramStruct,fourierStruct,K,n);
        ts = tinv([0.025  0.975],paramStruct.nblocks-1);      % 95 T-Score
        CI_c = simulated_call+ ts*scMC;
        CI_p = simulated_put+ ts*spMC;

        % Call option validation
        fprintf('nsteps = %3d, MC call price is %4.6f, std is %4.6f, abs err is %4.6f  ',n*25 ,simulated_call,scMC,abs(call_price-simulated_call));
        if (call_price > CI_c(1) && call_price < CI_c(2))
            fprintf('Call option is validated\n');
        else
            fprintf('Call option is not validated\n');
        end
        % Put option validation
        fprintf('nsteps = %3d, MC put price is %4.6f, std is %4.6f, abs err is %4.6f  ',n*25 ,simulated_put,spMC,abs(put_price-simulated_put));
        if (put_price > CI_p(1) && put_price < CI_p(2))
            fprintf(' Put option is validated\n');
        else
            fprintf('Put option is not validated\n');
        end
    end
end
fclose(Fid);
