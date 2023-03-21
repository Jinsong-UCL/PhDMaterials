function parBH = paramsWMSVtoBiHeston(t,V0,M,Q,R,beta,nFactors)
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Calibration and advanced simulation schemes for the Wishart Stochastic Volatility model
% (c) G. La Bua, D. Marazzina
% Corresponding Author: daniele.marazzina@polimi.it
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
A = [-M Q.'*Q; zeros(nFactors) M.'];

expA = expm(t*A);
expM_tran = expA(nFactors+1:2*nFactors, nFactors+1:2*nFactors);
expM = expM_tran.';

qt = expM * expA(1:nFactors, nFactors+1:2*nFactors);
muSigma = expM * V0 * expM_tran; % non-centrality matrix
%% eigenvaluee decomposition of qt
lambda1 = 0.5*(trace(qt)+sqrt(trace(qt)^2 - 4*det(qt)));
lambda2 = 0.5*(trace(qt)-sqrt(trace(qt)^2 - 4*det(qt)));

if qt(1,2)== 0 && qt(2,1) == 0
    X = eye(nFactors);
else
    v1 = [lambda1 - qt(2,2); qt(1,2)];
    v1 = v1/sqrt(v1.'*v1);
    v2 = [lambda2 - qt(2,2); qt(1,2)];
    v2 = v2/sqrt(v2.'*v2);
    
    X = [v1 v2];
end
tildeMu_X = X.'*muSigma*X;

zeta = diag(tildeMu_X);
epsilon = [lambda1;lambda2];

%% WSVM-BiHeston parameters mapping
v0Vec = diag(X.'*V0*X);
kappaVec = -log( zeta./ v0Vec) / t;
sigmaVec = 2*sqrt( epsilon .* kappaVec ./ (1-exp(-kappaVec * t)));
thetaVec = beta * sigmaVec.* sigmaVec ./ (4*kappaVec);

rhoVec = (diag(V0*R*Q)./(diag(V0).*sqrt(diag(Q.'*Q))));

parBH = [v0Vec(1),kappaVec(1),thetaVec(1),sigmaVec(1),rhoVec(1),v0Vec(2),kappaVec(2),thetaVec(2),sigmaVec(2),rhoVec(2)];

%EoF
end