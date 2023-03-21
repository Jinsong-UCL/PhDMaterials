function [v0,kappa,theta,eta,rho] = paramsWMSVtoHeston(t,V0,M,Q,R,beta,nFactors)
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Calibration and advanced simulation schemes for the Wishart Stochastic Volatility model
% (c) G. La Bua, D. Marazzina
% Corresponding Author: daniele.marazzina@polimi.it
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%>> Wishart-Heston mapping
% We assume that M is symmetric and negative definite

% initialize variables
A = [-M Q.'*Q; zeros(nFactors) M.'];

expA = expm(t*A);
expM_tran = expA(nFactors+1:2*nFactors, nFactors+1:2*nFactors);
expM = expM_tran.';

qt = expM * expA(1:nFactors, nFactors+1:2*nFactors);

muSigma = expM * V0 * expM_tran; % non-centrality matrix

% conditional moments of V_T and of Tr[V_T]
E1 = muSigma + beta*qt;
muTraceV = trace(E1);
varTraceV = 2*trace(beta*qt*qt + 2*muSigma*qt);

% Heston parameters
v0 = trace(V0);
if sum(sum(M-M.')) == 0 % M is symmetric
    theta = -0.5*beta*trace(M\Q.'*Q);
else
    theta = trace(sylv(M,M.',-beta*Q.'*Q));
end

kappa = abs(-log( (muTraceV-theta)/(trace(V0)-theta)))/t;
eta = abs(sqrt(varTraceV / ( (trace(V0)*exp(-kappa*t) + 0.5*theta*(1-exp(-kappa*t))) * (1-exp(-kappa*t))/kappa)));
rho = (trace(R*Q*V0) / sqrt(trace(V0)*trace(Q.'*Q*V0)));

%EoF
end