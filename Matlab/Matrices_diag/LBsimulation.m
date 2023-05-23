function [simulated_call, simulated_put,scMC,spMC, CF_e] = LBsimulation(market,param,fourier,K,n)
%% Retrieve parameters 
S0 = market.S0;
d = market.d;
T = market.T;

beta = param.beta;
An = param.An;
Am = param.Am;
Rn = param.Rn;
Rm = param.Rm;
R = Rn - Rm;
V_0 = param.V_0;

rho = param.rho;
kappa = param.kappa;
sigma = param.sigma;
sum(eig(sigma)>=0)

ngrid = fourier.ngrid; % number of grid points
xi = fourier.xi; 
%% Simulation parameters
% Number of simulations
nblocks = param.nblocks;
npaths =  param.npaths;
% Number of steps
nsteps = n*25;
dt = T/nsteps;

tic;