%% Warning option
[msgStr,warnId] = lastwarn;
warnStruct = warning('off',warnId);

%% Market parameters 
marketStruct.S0 = 1;
marketStruct.d = 2;
K_g = [50,75,90,100,110,125,150]/100;
T_g = [30,60,90,120,150,180,270,360]./360;
[K, T] = meshgrid(K_g,T_g);
% The following numbers are from Gnoatto and Grasselli 2014
%% Symmetric parameters
paramStruct.beta = 5.144;
paramStruct.An = [-0.0764 0.4837;0.4837 0.0739];
paramStruct.Am = [0.0679 0.6277;0.6277 -0.0520];
paramStruct.Rn = [-0.0725 0.0804;0.0804 0.0726];
paramStruct.Rm = [0.0941 0.0155;0.0155 -0.0941];
paramStruct.V_0 = [0.1688 0.0708;0.0708 0.1669];
% paramStruct.hn = -0.2218;
% paramStruct.hm = -0.1862;
paramStruct.hn = 0.05;
paramStruct.hm = 0.03;
paramStruct.rho = [-0.8417 0.1899;0.1170 -0.7834];
paramStruct.kappa = [2.2426,0.6764;0.0880,2.0778]; 
%paramStruct.sigma = [0.4368,0.1914;0.1914,0.7362]; %symmetric
paramStruct.sigma = [1.9368,0.2914;0.4966,0.5362]; %asymmetric
% Number of simulations
paramStruct.nblocks = 100;
paramStruct.npaths =  100;
%% Fourier parameters
fourierStruct.xwidth = 20; % width of the support in real space
fourierStruct.ngrid = 2^8; % number of grid points

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
 
%% Price grid
rate_1 = max(trace(paramStruct.Rn * paramStruct.V_0)+paramStruct.hn,0);
rate_2 = max(trace(paramStruct.Rm * paramStruct.V_0)+paramStruct.hm,0);

Call_price_grid_FT = zeros(length(T_g),length(K_g));
Put_price_grid_FT = zeros(length(T_g),length(K_g));
Call_price_grid_MC = zeros(length(T_g),length(K_g));
Put_price_grid_MC = zeros(length(T_g),length(K_g));
Impvol_FT =  zeros(length(T_g),length(K_g));
Impvol_MC =  zeros(length(T_g),length(K_g));

for t = 1:length(T_g)
    parfor k = 1:length(K_g)
        [msgStr,warnId] = lastwarn;
        warnStruct = warning('off',warnId);
        % Fourier Transform
        Call_price_grid_FT(t,k) = GCF(marketStruct,paramStruct,fourierStruct,K(k),T(t),1);
        Impvol_FT(t,k) = blsimpv(marketStruct.S0,K(k),rate_1,T(t),Call_price_grid_FT(t,k), 'Yield', rate_2,'Class', {'Call'},'Method','jackel2016');
        Put_price_grid_FT(t,k) = GCF(marketStruct,paramStruct,fourierStruct,K(k),T(t),-1);
        % Monte Carlo
        [Call_price_grid_MC(t,k), Put_price_grid_MC(t,k),~,~,~] = GGsimulation(marketStruct,paramStruct,fourierStruct,K(k),T(t),2);
        Impvol_MC(t,k) = blsimpv(marketStruct.S0,K(k),rate_1,T(t),Call_price_grid_MC(t,k), 'Yield', rate_2,'Class', {'Call'},'Method','jackel2016');
    end 
end
figure(1)
hold off
mesh(K,T,Call_price_grid_FT)
hold on
mesh(K,T,Call_price_grid_MC)
title('Call price surface for the matrix Heston model')
xlabel('K')
ylabel('T')
zlabel('Price')
legend('Fourier Transform','Monte Carlo')
savefig("Call price surface.fig")


figure(2)
hold off
mesh(K,T,Put_price_grid_FT)
hold on
surf(K,T,Put_price_grid_MC)
title('Put price surface for the matrix Heston model')
xlabel('K')
ylabel('T')
zlabel('Price')
legend('Fourier Transform','Monte Carlo')
savefig("Put price surface.fig")


figure(3)
hold off
mesh(K,T,Impvol_FT)
hold on
mesh(K,T,Impvol_MC)
title('Implied volatility surface for the matrix Heston model')
xlabel('K')
ylabel('T')
zlabel('Implied volatility')
legend('Fourier Transform','Monte Carlo')
savefig("Implied volatility.fig")

% % Call option
% figure(1)
% plot(xi,real(CF_E),xi,real(GCF_A),"--")
% axis([-20 20 -0.1 0.7])
% title('Real part of the discounted characteristic function for the matrix Heston model')
% xlabel('\xi')
% legend('Empirical','Analytical')
% savefig("Real part of the discounted characteristic function for the matrix Heston model(asymmetric).fig")
% 
% % Put option
% figure(2)
% plot(xi,imag(CF_E),xi,imag(GCF_A),"--")
% axis([-20 20 -0.3 0.3])
% title('Imaginary part of the discounted characteristic function for the matrix Heston model')
% xlabel('\xi')
% legend('Empirical','Analytical')
% savefig("Imaginary part of the discounted characteristic function for the matrix Heston model(asymmetric).fig")
