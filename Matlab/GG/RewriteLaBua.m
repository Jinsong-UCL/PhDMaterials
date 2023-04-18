S0 = 1;
K = 1;
% The following numbers are from Gnoatto and Grasselli 2014
am = [0.7764 0.4837;0.4837 0.9639];
an = [0.6679 0.6277;0.6277 0.8520];
hm = [0.2725 0.0804;0.0804 0.4726];
hn = [0.1841 0.0155;0.0155 0.4761];
y_0 = [0.1688 0.1708;0.1708 0.3169];
%am = [0.7764 0.0;0.0 0.9639];
%an = [0.6679 0.0;0.0 0.8520];
%hm = [0.2725 0.0;0.0 0.4726];
%hn = [0.1841 0.0;0.0 0.4761];
%y_0 = [0.1688 0.0;0.0 0.3169];
%y_bar =
rho = [-0.5417 0.1899;-0.1170 -0.4834];
kappa = [1.0426,0.6764;0.9880,0.8778]; 
sigma = [0.4364,0.1914;0.4966,0.7362];

M = -2*kappa;
Q = 2*sigma;
R = rho;
r = 0.02;
%rho = [-0.5417 0.1899;0.1899 -0.4834];
%kappa = [1.0426,0.6764;0.6764,0.8778]; 
%sigma = [0.4364,0.1914;0.1914,0.7362];

%rho = [-0.5417 0.0;0.0 -0.4834];
%kappa = [1.0426,0.0;0.0,0.8778]; 
%sigma = [0.4364,0.0;0.0,0.7362];

h = hm - hn;
beta = 3.1442;
T = 1;
d = 2;

nsteps = 100;
nsims = 100;
dt = T/nsteps;
sqrt_dt = sqrt(dt);

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
yPaths = y_0(:,:,:,ones(1,nsims));

%% Data for (log-) asset simulation

%%% Initialization
XPaths = zeros(nsteps+1,nsims);
XPaths(1,:) = log(S0);

%% j-loop (time), i-loop (Simulations) and l-loop (Assets)
%Vt_mean = zeros(2);
for step = 1:nsteps
    
    % Generate random numbers
    NormRV = randn(nsims, 7); % normal random variables
    NormRVdt = NormRV(:,5:7) * sqrt_dt; % normal random variables multiplied by square root of delta t
    UnifRV = rand(nsims, 2); % uniform random variables
    GammaRV = gamrnd_fast((beta-1)/2, 2, nsims, 2); % gamma random variables
    
    for sim = 1:nsims
        Vt = yPaths(:,:,:, sim); % get the current V path for this simulation
        
        %%% Exact Wishart process simulation (Ahdida-Alfonsi) P1034 Eq.21 
        modX = f1*Vt*f2; % calculate modX based on the current V path
        u11 = modX(1,1) - modX(1,2)^2 / modX(2,2); % calculate u11
        u12 = modX(1,2) / sqrt(modX(2,2)); % calculate u12
        u22 = modX(2,2); % calculate u22
    
        % calculate U11, U12, and U22 based on modX
        lambda = u11/dt;
        Z = lambda + 2*log(UnifRV(sim,1));
        if Z > 0
            U11 = (GammaRV(sim,1) + ((NormRV(sim,1) + sqrt(Z))^2 + NormRV(sim,2)^2))*dt;
        else
            U11 = GammaRV(sim,1)*dt;
        end
        U12 = u12 + NormRVdt(sim,1);
        U22 = u22;
        Y = [U11+U12^2 U12*sqrt(U22); U12*sqrt(U22) U22];
        
        % calculate Y based on P and U11, U12, and U22
        Y = P*Y*P;
        u11 = Y(1,1) - Y(1,2)^2 / Y(2,2); % calculate new u11
        u12 = Y(1,2) / sqrt(Y(2,2)); % calculate new u12
        u22 = Y(2,2); % calculate new u22
        
        % calculate U11, U12, and U22 based on new u11, u12, and u22
        lambda = u11/dt;
        Z = lambda + 2*log(UnifRV(sim,2));
        if Z > 0
            U11 = (GammaRV(sim,2) + ((NormRV(sim,3) + sqrt(Z))^2 + NormRV(sim,4)^2))*dt;
        else
            U11 = GammaRV(sim,2)*dt;
        end
        U12 = u12 + NormRVdt(sim,2);
        U22 = u22;
        Y = [U11+U12^2 U12*sqrt(U22); U12*sqrt(U22) U22];
        
        % calculate new Y based on P and U11, U12, and U22
        Y = P*Y*P;
        Vt_dt = theta*Y*thetat; % calculate the new Vt_dt for this simulation
        
        yPaths(:,:,:, sim) = Vt_dt; % set the new V path for this simulation
        
        %% Asset Simulation
        % Calculate intermediate variables for asset simulation
        traceVt = trace(Vt);
        traceVt_dt = trace(Vt_dt);
        dTrace_dt = traceVt_dt - traceVt;
        rho_t = trace(R*Q*Vt) / (sqrt(traceVt)*sqrt(trace(Q.'*Q*Vt)));
        
        % Limit rho_t to be within [-1,1]
        rho_t = min(max(rho_t,-1),1);
        
        % Calculate asset path
        w_t = 0.5* (dTrace_dt - ( beta*trace(Q.'*Q) + 2*trace(M*Vt))*dt) / sqrt(trace(Vt*Q.'*Q));
      
        XPaths(step+1,sim) = XPaths(step,sim) - 0.5*traceVt*dt + rho_t*sqrt(traceVt)*w_t + sqrt((1-rho_t^2)*traceVt)*NormRVdt(sim,3);
        
    end
    
end

SPaths = exp((r)*T) * exp(XPaths);

function r = gamrnd_fast(a, b, nRow, nColumn)
    r = b .* randg(a,[nRow, nColumn]);
end