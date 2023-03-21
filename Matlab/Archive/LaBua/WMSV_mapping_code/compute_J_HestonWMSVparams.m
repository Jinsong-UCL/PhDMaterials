function J = compute_J_HestonWMSVparams(t,V0,M,Q,R,beta,nFactors)
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

% conditional moments of V_T and of Tr[V_T]
E1 = trace(muSigma + beta*qt);
muTraceV = trace(E1);
varTraceV = 2*trace(beta*qt*qt + 2*muSigma*qt);

% Heston parameters
v0 = trace(V0);
theta = -0.5*beta*trace(M\Q.'*Q);
kappa = -log( (muTraceV-theta)/(trace(V0)-theta))/t;
sigma = sqrt(varTraceV / ( (trace(V0)*exp(-kappa*t) + 0.5*theta*(1-exp(-kappa*t))) * (1-exp(-kappa*t))/kappa));

denomSigma = ( (trace(V0)*exp(-kappa*t) + 0.5*theta*(1-exp(-kappa*t))) * (1-exp(-kappa*t))/kappa);
% Initialize J
nParams = (nFactors*(nFactors+1)) + 2*nFactors^2 + 1; % V0 and M symmetric
J = zeros(nParams,5);

%% loop for V0 (symmetric matrix)
for i = 1:nFactors
    
    for j = i:nFactors 
        %%% derivative of mu_Sigma wrt V(i,j)
        auxV0 = zeros(nFactors);
        auxV0(i,j) = 1;
        dMu_Vij = expM * auxV0 * expM_tran;
        
        kron = 1; %kronecker's delta for i = j
        if i ~= j
            dMu_Vij = dMu_Vij + expM * auxV0.' * expM_tran;
            kron = 0;
        end
        
        %%% derivative of kHeston wrt V(i,j)
        dkappaH_Vij = -1/t * ( (trace(V0) - theta)*trace(dMu_Vij)  -  (E1-theta) * kron)  / ( (E1-theta)*(trace(V0)-theta));
        
        %%% derivative of sigmaHeston wrt V(i,j)
        
        % derivative of quantity that appears in the denominator of sigma
        dDen_Vij = (t*exp(-2*kappa*t)/kappa * (2*trace(V0)-theta - exp(t*kappa)*(trace(V0)-theta)) - denomSigma / kappa) * dkappaH_Vij - exp(-2*kappa*t)*(1-exp(kappa*t))/kappa*kron;
        
        dsigmaH_Vij = 0.5/sigma  * ((4*trace(dMu_Vij*qt)*denomSigma - varTraceV * dDen_Vij) / denomSigma^2);
        
        %%% derivative of rhoHeston wrt V(i,j)
        drhoH_Vij = trace(R*Q*auxV0)/sqrt(trace(V0)*trace(Q.'*Q*V0)) - 0.5*trace(R*Q*V0)*(trace(Q.'*Q*V0)*kron + trace(V0)*trace(Q.'*Q*auxV0)) / (trace(V0)*trace(Q.'*Q*V0))^(3/2);
        
        if i ~= j
            drhoH_Vij = drhoH_Vij +  trace(R*Q*auxV0.')/sqrt(trace(V0)*trace(Q.'*Q*V0)) - 0.5*trace(R*Q*V0)* trace(V0)*trace(Q.'*Q*auxV0.') / (trace(V0)*trace(Q.'*Q*V0))^(3/2); % giusto
        end
        
        %%% allocate derivative wrt V0(i,j) in matrix J
        index = (i-1)*(nFactors - 0.5*i) + j; % mapping from upper triangular matrix to vector
        
        if j == i
            % d v0(i) / d V(i,i)
            J(index,1) = 1;
        end
        
        J(index,2) = dkappaH_Vij;
        J(index,4) = dsigmaH_Vij;
        J(index,5) = drhoH_Vij;
        
    end
    
end

%% loop for M (non-symmetric matrix)
startM = nFactors * (nFactors + 1)/2; %position after which the elements in M are placed in matrix J
for i = 1:nFactors
    for j = i:nFactors 
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
        
        aux1M = zeros(nFactors);
        aux1M(i,j) = 1;
        
        %%% derivative of thetaHeston wrt M(i,j)
        dthetaH_Mij = 0.5*beta*trace(M\ aux1M *(M\Q.'*Q));
        
        if i~=j
            dthetaH_Mij = 2*dthetaH_Mij;
        end
        
        %%% derivative of kHeston wrt M(i,j)
        dkappaH_Mij = -1/t  * ( (trace(V0) - theta)* trace(dMu_Mij + beta*dqt_Mij) + (E1- trace(V0)) * dthetaH_Mij)  / ...
            ( (E1-theta)*(trace(V0)-theta));
        
        %%% derivative of sigmaHeston wrt M(i,j)
        
        % derivative of quantity that appears in the denominator of sigma
        dDen_Mij = (t*exp(-2*kappa*t)/kappa * (2*trace(V0)-theta - exp(t*kappa)*(trace(V0)-theta)) - denomSigma / kappa) * dkappaH_Mij +...
            (1-exp(-kappa*t))^2 /(2*kappa)  * dthetaH_Mij;
        dsigmaH_Mij = 0.5/sigma  * ((4*trace((muSigma+beta*qt)*dqt_Mij + qt * dMu_Mij )*denomSigma - varTraceV * dDen_Mij) / denomSigma^2);
        
        %%% allocate derivative wrt M(i,j) in matrix J
        index = (i-1)*(nFactors - 0.5*i) + j; % mapping from upper triangular matrix to vector
        index = startM + index;
        
        J(index,2) = dkappaH_Mij;
        J(index,3) = dthetaH_Mij;
        J(index,4) = dsigmaH_Mij;
        
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
        
        %%% derivative of thetaHeston wrt Q(i,j)
        dthetaH_Qij = -0.5*beta * trace(M\ (aux1Q*Q + (aux1Q*Q).'));
        
        %%% derivative of kHeston wrt Q(i,j)
        dkappaH_Qij = -1/t  *  (beta * (trace(V0) - theta)*trace(dqt_Qij) +...
            (E1- trace(V0)) * dthetaH_Qij )  / ( (E1-theta)*(trace(V0)-theta));
        
        %%% derivative of sigmaHeston wrt Q(i,j)
        
        % derivative of quantity that appears in the denominator of sigma
        dDen_Qij = (t*exp(-2*kappa*t)/kappa * (2*trace(V0)-theta - exp(t*kappa)*(trace(V0)-theta)) - denomSigma / kappa) * dkappaH_Qij +...
            (1-exp(-kappa*t))^2 /(2*kappa)  * dthetaH_Qij;
        dsigmaH_Qij = 0.5/sigma  * ((4*trace((muSigma+beta*qt)*dqt_Qij)*denomSigma - varTraceV * dDen_Qij) / denomSigma^2);
        
        %%% derivative of rhoHeston wrt Q(i,j)
        drhoH_Qij = trace(R*aux1Q.'*V0) / sqrt(trace(V0) * trace(Q.'*Q*V0)) -...
            trace(R*Q*V0)*trace(Q.'*aux1Q.' * V0) / (sqrt(trace(V0)) * trace(Q.'*Q*V0)^(3/2));
        
        %%% allocate derivative wrt Q(i,j) in matrix J
        
        index = startQ + j + (i-1)*nFactors; % mapping from matrix Q to J
        
        J(index,2) = dkappaH_Qij;
        J(index,3) = dthetaH_Qij;
        J(index,4) = dsigmaH_Qij;
        J(index,5) = drhoH_Qij;
    end
end
%% loop for R (non-symmetric matrix)
startR = startQ + nFactors^2; %position after which the elements in R are placed in matrix J

for i = 1:nFactors
    for j = 1:nFactors
        auxR = zeros(nFactors);
        auxR(i,j) = 1;
        
        drhoH_Rij = trace(auxR*Q*V0) / sqrt(trace(V0)*trace(Q.'*Q*V0));
        
        %%% allocate derivative wrt R(i,j) in matrix J
        index = startR + j + (i-1)*nFactors; % mapping from matrix R to J
        
        J(index,5) = drhoH_Rij;
    end
end

%% derivatives wrt to beta

%%% derivative of thetaHeston wrt beta)
dthetaH_beta = -0.5*trace(M\Q.'*Q);

%%% derivative of kHeston wrt beta
dkappaH_beta = -1/t  * (trace(qt) * (trace(V0) - theta) + (E1-trace(V0)) * dthetaH_beta) / ( (E1-theta)*(trace(V0)-theta));

%%% derivative of sigmaHeston wrt M(i,j)

% derivative of quantity that appears in the denominator of sigma
dDen_beta = (t*exp(-2*kappa*t)/kappa * (2*trace(V0)-theta - exp(t*kappa)*(trace(V0)-theta)) - denomSigma / kappa) * dkappaH_beta + (1-exp(-kappa*t))^2 /(2*kappa)  * dthetaH_beta;
dsigmaH_beta = 0.5/sigma  * ((2*trace(qt*qt)*denomSigma - varTraceV * dDen_beta) / denomSigma^2);

%%% allocate derivative wrt beta in matrix J
J(end,2) = dkappaH_beta;
J(end,3) = dthetaH_beta;
J(end,4) = dsigmaH_beta;

%EoF
end