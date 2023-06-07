%% Warning option
[msgStr,warnId] = lastwarn;
warnStruct = warning('off',warnId);

%% Market parameters 
marketStruct.S0 = 100;
marketStruct.d = 2;
K_g = 50:10:150;
T_g = 0:0.1:1;
[K, T] = meshgrid(K_g,T_g);
% The following numbers are from Gnoatto and Grasselli 2014
%% Symmetric parameters
paramStruct.beta = 3.144;
paramStruct.An = [0.7764 0.4837;0.4837 0.9639];
paramStruct.Am = [0.6679 0.6277;0.6277 0.8520];
paramStruct.Rn = [0.2725 0.0804;0.0804 0.4726];
paramStruct.Rm = [0.1841 0.0155;0.0155 0.4761];
paramStruct.V_0 = [0.1688 0.1708;0.1708 0.3169];
paramStruct.hn = -0.2218;
paramStruct.hm = -0.1862;
paramStruct.rho = [-0.5417 0.1899;-0.1170 -0.4834];
paramStruct.kappa = [1.0426,0.6764;0.9880,0.8778]; 
%paramStruct.sigma = [0.4368,0.1914;0.1914,0.7362]; %symmetric
paramStruct.sigma = [0.4368,0.1914;0.4966,0.7362]; %asymmetric
% Number of simulations
paramStruct.nblocks = 10;
paramStruct.npaths =  10;
%% Fourier parameters
fourierStruct.xwidth = 20; % width of the support in real space
fourierStruct.ngrid = 2^7; % number of grid points

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
rate = max(trace(paramStruct.Rn * paramStruct.V_0)+paramStruct.hn,0);

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
        Impvol_FT(t,k) = blsimpv(marketStruct.S0,K(k),rate,T(t),Call_price_grid_FT(t,k));
        Put_price_grid_FT(t,k) = GCF(marketStruct,paramStruct,fourierStruct,K(k),T(t),-1);
        % Monte Carlo
        [Call_price_grid_MC(t,k), Put_price_grid_MC(t,k),~,~,~] = GGsimulation(marketStruct,paramStruct,fourierStruct,K(k),T(t),2);
        Impvol_MC(t,k) = blsimpv(marketStruct.S0,K(k),rate,T(t),Call_price_grid_MC(t,k));
    end 
end
figure(1)
hold off
mesh(K,T,Call_price_grid_FT)
%hold on
%surf(K,T,Call_price_grid_MC,'FaceColor','b')
title('Call price surface for the matrix Heston model')
xlabel('K')
ylabel('T')
zlabel('Price')
legend('Fourier Transform','Monte Carlo')
savefig("Call price surface.fig")


figure(2)
surf(K,T,Put_price_grid_FT,'FaceColor','r')
%hold on
%surf(K,T,Put_price_grid_MC,'FaceColor','b')
title('Put price surface for the matrix Heston model')
xlabel('K')
ylabel('T')
zlabel('Price')
legend('Fourier Transform','Monte Carlo')
savefig("Put price surface.fig")


figure(3)
mesh(K,T,Impvol_FT,'FaceColor','r')
hold on
mesh(K,T,Impvol_MC,'FaceColor','b')
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
