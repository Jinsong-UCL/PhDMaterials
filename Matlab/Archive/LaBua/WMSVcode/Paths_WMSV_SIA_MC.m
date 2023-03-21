function SPaths = Paths_WMSV_SIA_MC(params, parH, S0, t, tSteps, nSims)
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Calibration and advanced simulation schemes for the Wishart Stochastic Volatility model
% (c) G. La Bua, D. Marazzina
% Corresponding Author: daniele.marazzina@polimi.it
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%>> SIA Monte Carlo Algorithm

V0 = params.Sigma0;
Q = params.Q;
M = params.M;
R = params.R;
beta = params.beta;

r = params.r;
q = params.q;

if beta < 1
    error('wishart:wishart3:betaNotValid', 'Beta is not valid');
end

dt = t/tSteps;
sqrt_dt = sqrt(dt);

%% parameters of approximating Heston displaced model
kappaH = parH(1);
thetaH = parH(2);
etaH = parH(3);

%% Data for Wishart Process (Variance-Covariance Matrix)
%%% integral of matrix exponential

ret = expm([-M Q'*Q; zeros(2) M']*dt);
fact1 = ret(3:4, 3:4);
fact2 = ret(1:2, 3:4);

qt = fact1.' * fact2;

mt = expm(dt*M);

theta = chol(qt/dt, 'lower');
thetat = theta.';
itheta = inv(theta).';
f1 = theta\mt;
f2 = mt.' * itheta;

%%% permutation Matrix
P = [0 1; 1 0];

%%% Initialization
VPaths = V0(:,:,:,ones(1,nSims));

%% Data for (log-) asset simulation

% Initialization
XPaths = zeros(tSteps+1,nSims);
XPaths(1,:) = log(S0);

%% j-loop (time), i-loop (Simulations) and l-loop (Assets)
%Vt_mean = zeros(2);
for j = 1:tSteps
     
    % random numbers generation
    NormRV = randn(nSims, 7);
    NormRVdt = NormRV(:,5:6) * sqrt_dt;
    UnifRV = rand(nSims, 2);
    GammaRV = gamrnd_fast((beta-1)/2, 2, nSims, 2);
    
    for i = 1:nSims
        Vt = VPaths(:,:,:, i);
        
        %% Exact Wishart process simulation (Ahdida-Alfonsi)     
        modX = f1*Vt*f2;
        u11 = modX(1,1) - modX(1,2)^2 / modX(2,2);
        u12 = modX(1,2) / sqrt(modX(2,2));
        u22 = modX(2,2);
        
        %
        lambda = u11/dt;
        Z = lambda + 2*log(UnifRV(i,1));
        if Z > 0
            U11 = (GammaRV(i,1) + ((NormRV(i,1) + sqrt(Z))^2 + NormRV(i,2)^2))*dt;
        else
            U11 = GammaRV(i,1)*dt;
        end
        U12 = u12 + NormRVdt(i,1);
        U22 = u22;
        Y = [U11+U12^2 U12*sqrt(U22); U12*sqrt(U22) U22];
        
        %
        Y = P*Y*P;
        u11 = Y(1,1) - Y(1,2)^2 / Y(2,2);
        u12 = Y(1,2) / sqrt(Y(2,2));
        u22 = Y(2,2);
        
        %
        lambda = u11/dt;
        Z = lambda + 2*log(UnifRV(i,2));
        if Z > 0
            U11 = (GammaRV(i,2) + ((NormRV(i,3) + sqrt(Z))^2 + NormRV(i,4)^2))*dt;
        else
            U11 = GammaRV(i,2)*dt;
        end
        U12 = u12 + NormRVdt(i,2);
        U22 = u22;
        Y = [U11+U12^2 U12*sqrt(U22); U12*sqrt(U22) U22];
        
        %
        Y = P*Y*P;
        Vt_dt = theta*Y*thetat;
        
        VPaths(:,:,:, i) = Vt_dt;
        
        %% Asset Simulation
        
        traceVt = trace(Vt);
        traceVt_dt = trace(Vt_dt);
                
        rho_t = trace(R*Q*Vt) / (sqrt(traceVt)*sqrt(trace(Q.'*Q*Vt)));
        
        rho_t = min(max(rho_t,-1),1);
        
        Z_x = NormRV(i,7);
        K0 = -dt*rho_t*kappaH*thetaH/etaH;
        % gamma1 e gamma2 = 0.5
        K1 = 0.5*dt*(kappaH*rho_t/etaH - 0.5) - rho_t / etaH;
        K2 = 0.5*dt*(kappaH * rho_t / etaH - 0.5) + rho_t / etaH;
        K3 = 0.5*dt*(1-rho_t^2);
        K4 = 0.5*dt*(1-rho_t^2); % potrei anche mettere rho_t_dt
        XPaths(j+1,i) = XPaths(j,i) + K0 + K1* traceVt + K2*traceVt_dt + sqrt(K3*traceVt + K4*traceVt_dt) * Z_x;
        
    end
    
end

SPaths = exp((r-q)*t) * exp(XPaths);
%EoF
end

function r = gamrnd_fast(a, b, nRow, nColumn)

r = b .* randg(a,[nRow, nColumn]);
%EoF
end
