function J = compute_J_BiHestonWMSVparams(t,V0,M,Q,R,beta,nFactors)
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Calibration and advanced simulation schemes for the Wishart Stochastic Volatility model
% (c) G. La Bua, D. Marazzina
% Corresponding Author: daniele.marazzina@polimi.it
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% We assume that M is symmetric and negative definite

% initialize variables
A = [-M Q.'*Q; zeros(nFactors) M.'];

expA = expm(t*A);
expM_tran = expA(nFactors+1:2*nFactors, nFactors+1:2*nFactors);
expM = expM_tran.';

qt = expM * expA(1:nFactors, nFactors+1:2*nFactors);

muSigma = expM * V0 * expM_tran; % non-centrality matrix

% Bi-Heston parameters
parBH = paramsWSVMtoBiHeston(t,V0,M,Q,R,beta,nFactors);

v0(:,1) = parBH([1,6]);
kappa(:,1) = parBH([2,7]);
theta(:,1) = parBH([3,8]);
sigma(:,1) = parBH([4,9]);

%% eigenvaluee decomposition of qt
lambda1 = 0.5*(trace(qt)+sqrt(trace(qt)^2 - 4*det(qt)));
lambda2 = 0.5*(trace(qt)-sqrt(trace(qt)^2 - 4*det(qt)));

if qt(1,2) == 0 && qt(2,1) == 0
    X = eye(nFactors);
    v1 = [1; zeros(nFactors-1,1)];
    norm1 = 1;
    
    v2 = [zeros(nFactors-1,1);1];
    norm2 = 1;
    
    DiagonalCheck = 1;
else
    v1 = [lambda1 - qt(2,2); qt(1,2)];
    norm1 = sqrt(v1.'*v1);
    v1 = v1/norm1;
    v2 = [lambda2 - qt(2,2); qt(1,2)];
    norm2 = sqrt(v2.'*v2);
    v2 = v2/norm2;
    
    X = [v1 v2];
    
    DiagonalCheck = 0;
end

tildeMu_X = X.'*muSigma*X;

zeta = diag(tildeMu_X);
epsilon = [lambda1;lambda2];

% useful quantities
RQ = R*Q;
QQ = Q.'*Q;
denRho = diag(V0).*sqrt(diag(Q.'*Q));

%% Initialize J
nParams = (nFactors*(nFactors+1)) + 2*nFactors^2 + 1; % V0 and M symmetric
J = zeros(nParams,10);

%% loop for V0 (symmetric matrix)
for i = 1:nFactors
    
    for j = i:nFactors % as symmetric matrix
        %%% derivative of mu_Sigma wrt V(i,j)
        auxV0 = zeros(nFactors);
        auxV0(i,j) = 1;
        dMu_Vij = expM * auxV0 * expM_tran;
        
        if i ~= j
            dMu_Vij = dMu_Vij + expM * auxV0.' * expM_tran;
        end
        
        %%% derivative of v0 wrt V(i,j)
        dv0BH_Vij = diag(X.'*auxV0*X);
        
        if i~= j
            dv0BH_Vij = dv0BH_Vij + diag(X.'*auxV0.'*X);
        end
        %%% derivative of kappaBiHeston wrt V(i,j)
        dzeta_Vij = diag(X.'*dMu_Vij*X);
        
        dkappaBH_Vij = 1/t * (dv0BH_Vij ./ v0 - dzeta_Vij ./zeta);
        
        %%% derivative of sigmaBiHeston wrt V(i,j)
        
        dsigmaBH_Vij = 0.5*sigma.*dkappaBH_Vij./kappa .* (1+t*kappa-exp(t*kappa))./(1-exp(t*kappa));
        
        %%% derivative of thetaBiHeston wrt V(i,j)
        
        dthetaBH_Vij = theta./ (sigma.*kappa) .* (2*kappa.* dsigmaBH_Vij - sigma.*dkappaBH_Vij);
        %%% derivative of rhoHeston wrt V(i,j)
        drhoBH_Vij = zeros(nFactors,1);
        
        drhoBH_Vij_temp = flipud(diag(fliplr(RQ)))./denRho;
        
        if i == j
            drhoBH_Vij(i) = -drhoBH_Vij_temp(i)*V0(1,2)/V0(i,i);
        else
            drhoBH_Vij = drhoBH_Vij_temp;
        end
        %%% allocate derivative wrt V0(i,j) in matrix J
        index = (i-1)*(nFactors - 0.5*i) + j; % mapping from upper triangular matrix to vector
        
        J(index,[1,6]) = dv0BH_Vij;
        J(index,[2,7]) = dkappaBH_Vij;
        J(index,[3,8]) = dthetaBH_Vij;
        J(index,[4,9]) = dsigmaBH_Vij;
        J(index,[5,10]) = drhoBH_Vij;
        
    end
    
end

%% loop for M (symmetric matrix)
startM = nFactors * (nFactors + 1)/2; %position after which the elements in M are placed in matrix J
for i = 1:nFactors
    for j = i:nFactors % as symmetric matrix
        dA_Mij = zeros(2*nFactors);
        
        dA_Mij(i,j) = -1;
        dA_Mij(nFactors+j,nFactors+i) = 1;
        
        C_M = [A dA_Mij; zeros(2*nFactors) A];
        D_M = expm(t*C_M);
        dF_Mij = D_M(1:2*nFactors,2*nFactors+1:4*nFactors);
        
        E1_Mij = dF_Mij(nFactors+1:2*nFactors,nFactors+1:2*nFactors);
        E2_Mij = dF_Mij(1:nFactors,nFactors+1:2*nFactors);
        
        %%% derivative of mu_Sigma wrt M(i,j)
        dMu_Mij = expM * V0 * E1_Mij + (expM * V0 * E1_Mij).';
        
        %%% derivative of qt wrt M(i,j)
        dqt_Mij = E1_Mij.' * expA(1:nFactors,nFactors+1:2*nFactors) + expM * E2_Mij;
        
        %%% adjust in order to take into account that M is symmetric
        if i~= j
            dA_Mij = dA_Mij.';
            
            C_M = [A dA_Mij; zeros(2*nFactors) A];
            D_M = expm(t*C_M);
            dF_Mij = D_M(1:2*nFactors,2*nFactors+1:4*nFactors);
            
            E1_Mij = dF_Mij(nFactors+1:2*nFactors,nFactors+1:2*nFactors);
            E2_Mij = dF_Mij(1:nFactors,nFactors+1:2*nFactors);
            
            dMu_Mij = dMu_Mij + expM * V0 * E1_Mij + (expM * V0 * E1_Mij).';
            dqt_Mij = dqt_Mij + E1_Mij.' * expA(1:nFactors,nFactors+1:2*nFactors) + expM * E2_Mij;
            
        end
        
        %%% derivative of matrix X wrt M(i,j)
        if DiagonalCheck == 1
            auxX = zeros(nFactors);
            dv1_Mij = zeros(nFactors,1);
            dv2_Mij = zeros(nFactors,1);
            
            dX_Mij = zeros(nFactors);
        else
            auxX = (trace(qt)*trace(dqt_Mij) - 2*det(qt)*trace(qt\dqt_Mij))/sqrt(trace(qt)^2 - 4*det(qt));
            dv1_Mij = [0.5*(trace(dqt_Mij)+auxX)-dqt_Mij(2,2);dqt_Mij(1,2)];
            dv2_Mij = [0.5*(trace(dqt_Mij)-auxX)-dqt_Mij(2,2);dqt_Mij(1,2)];
            dX_Mij = [(dv1_Mij - v1.*(v1.'*dv1_Mij))/norm1,  (dv2_Mij - v2.*(v2.'*dv2_Mij))/norm2];
        end
        
        %%% derivative of v0BiHeston wrt M(i,j)
        dv0BH_Mij = diag(dX_Mij.'*V0*X + X.'*V0*dX_Mij);
        
        %%% derivative of kBiHeston wrt M(i,j)
        dzeta_Mij = diag(dX_Mij.'*muSigma*X + X.'*dMu_Mij*X + X.'*muSigma*dX_Mij);
        
        dkappaBH_Mij = (dv0BH_Mij./v0 - dzeta_Mij./zeta)/t;
        
        %%% derivative of sigmaBiHeston wrt M(i,j)
        depsilon_Mij = [dv1_Mij(1)+dqt_Mij(2,2);dv2_Mij(1)+dqt_Mij(2,2)];
        
        dsigmaBH_Mij = 0.5*sigma.*(depsilon_Mij./epsilon + (1./kappa + t./(1-exp(kappa*t))).*dkappaBH_Mij);
        
        %%% derivative of thetaBiHeston wrt M(i,j)
        dthetaBH_Mij = theta.*(2*kappa.*dsigmaBH_Mij - sigma .* dkappaBH_Mij)./ (kappa.*sigma);
        
        %%% allocate derivative wrt M(i,j) in matrix J
        index = (i-1)*(nFactors - 0.5*i) + j; % mapping from upper triangular matrix to vector
        index = startM + index;
        
        J(index,[1,6]) = dv0BH_Mij;
        J(index,[2,7]) = dkappaBH_Mij;
        J(index,[3,8]) = dthetaBH_Mij;
        J(index,[4,9]) = dsigmaBH_Mij;
    end
end

%% loop for Q (non-symmetric matrix)
startQ = startM + nFactors * (nFactors + 1)/2; %position after which the elements in Q are placed in matrix J (considered that M is symmetric!)
for i = 1:nFactors
    for j = 1:nFactors
        
        dA_Qij = zeros(2*nFactors);
        
        aux1Q = zeros(nFactors);
        aux1Q(j,i) = 1;
        dA_Qij(1:nFactors,nFactors+1:2*nFactors) = aux1Q*Q + (aux1Q*Q).';
        
        C_Q = [A dA_Qij; zeros(2*nFactors) A];
        D_Q = expm(t*C_Q);
        dF_Qij = D_Q(1:2*nFactors,2*nFactors+1:4*nFactors);
        
        %E1_Qij = dF_Qij(nFactors+1:2*nFactors,nFactors+1:2*nFactors);
        E2_Qij = dF_Qij(1:nFactors,nFactors+1:2*nFactors);
        
        %%% derivative of qt wrt Q(i,j)
        dqt_Qij = expM * E2_Qij;
        
        %%% derivative of matrix X wrt Q(i,j)
        if DiagonalCheck == 1
            auxX = zeros(nFactors);
            
            dv1_Qij = zeros(nFactors,1);
            dv2_Qij = zeros(nFactors,1);
            
            dX_Qij = zeros(nFactors);
            
        else
            auxX = (trace(qt)*trace(dqt_Qij) - 2*det(qt)*trace(qt\dqt_Qij))/sqrt(trace(qt)^2 - 4*det(qt));
            
            dv1_Qij = [0.5*(trace(dqt_Qij)+auxX)-dqt_Qij(2,2);dqt_Qij(1,2)];
            dv2_Qij = [0.5*(trace(dqt_Qij)-auxX)-dqt_Qij(2,2);dqt_Qij(1,2)];
            
            dX_Qij = [(dv1_Qij - v1.*(v1.'*dv1_Qij))/norm1,  (dv2_Qij - v2.*(v2.'*dv2_Qij))/norm2];
            
        end
        
        
        
        %%% derivative of v0BiHeston wrt Q(i,j)
        dv0BH_Qij = diag(dX_Qij.'*V0*X + X.'*V0*dX_Qij);
        
        %%% derivative of kBiHeston wrt Q(i,j)
        dzeta_Qij = diag(dX_Qij.'*muSigma*X + X.'*muSigma*dX_Qij);
        
        dkappaBH_Qij = (dv0BH_Qij./v0 - dzeta_Qij./zeta)/t;
        
        %%% derivative of sigmaBiHeston wrt Q(i,j)
        depsilon_Qij = [dv1_Qij(1)+dqt_Qij(2,2);dv2_Qij(1)+dqt_Qij(2,2)];
        
        dsigmaBH_Qij = 0.5*sigma.*(depsilon_Qij./epsilon + (1./kappa + t./(1-exp(kappa*t))).*dkappaBH_Qij);
        
        %%% derivative of thetaBiHeston wrt Q(i,j)
        dthetaBH_Qij = theta.*(2*kappa.*dsigmaBH_Qij - sigma .* dkappaBH_Qij)./ (kappa.*sigma);        %%% allocate derivative wrt Q(i,j) in matrix J
        
        %%% derivative of rhoBiHeston wrt Q(i,j)
        VRe = V0*R*aux1Q.';
        VRQ = V0*R*Q;
        
        drhoBH_Qij = (diag(VRe).*diag(V0).*sqrt(diag(QQ)) - 0.5*diag(VRQ).*diag(V0).*diag(aux1Q*Q+Q.'*aux1Q.')./diag(sqrt(QQ)))./ (diag(V0).^2 .* diag(QQ));
        
        %%% allocate derivative wrt M(i,j) in matrix J
        index = startQ + j + (i-1)*nFactors; % mapping from matrix Q to J
        
        J(index,[1,6]) = dv0BH_Qij;
        J(index,[2,7]) = dkappaBH_Qij;
        J(index,[3,8]) = dthetaBH_Qij;
        J(index,[4,9]) = dsigmaBH_Qij;
        J(index,[5,10]) = drhoBH_Qij;
    end
end
%% loop for R (non-symmetric matrix)
startR = startQ + nFactors^2; %position after which the elements in R are placed in matrix J

for i = 1:nFactors
    for j = 1:nFactors
        
        drhoBH_Rij = Q(j,:).'.*V0(:,i)./(sqrt(diag(QQ)).*diag(V0));
        
        %%% allocate derivative wrt R(i,j) in matrix J
        index = startR + j + (i-1)*nFactors; % mapping from matrix R to J
        
        J(index,[5,10]) = drhoBH_Rij;
    end
end

%% derivatives wrt to beta

%%% derivative of thetaHeston wrt beta)
dthetaBH_beta = theta/beta;

%%% allocate derivative wrt beta in matrix J
J(end,[3,8]) = dthetaBH_beta;

%EoF
end