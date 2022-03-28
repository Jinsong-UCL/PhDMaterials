% Set process parameters, Schoutens page 82
function parameters = parameters(distr,T,dt,rf,q)

parameters.distr = distr;
parameters.T = T;
parameters.dt = dt;
parameters.rf = rf;
if nargin == 5
    parameters.q = q;
else
    parameters.q = 0;
end

switch distr

    case 1 % Normal

        m = 0; % mean
        %s = 0.17801; % standard deviation
        s = 0.4; % standard deviation
        %s = 0.3; % standard deviation

        % Rearrange parameters (time rescaling)
        parameters.m = m;
        parameters.s = s*sqrt(dt);

        % Compute vol(T), Schoutens page 60
        %parameters.vol = s*sqrt(T);
        
        % analiticity bounds
        % 0 if \pm\infty
        parameters.lambdam = 0;
        parameters.lambdap = 0;

    case 2 % NIG, Schoutens page 59

        %alpha = 6.1882;
        %beta = -3.8941;
        %delta = 0.1622;
        alpha = 15;
        beta = -5;
        delta = 0.5;

        % Rearrange parameters (time rescaling)
        parameters.alpha = alpha;
        parameters.beta = beta;
        parameters.delta = delta*dt;

        % Compute vol(T), Schoutens page 60
        c = alpha^2 - beta^2;
        variance = alpha^2*delta*c^(-1.5);
        %parameters.vol = sqrt(variance*T);
        
        % analiticity bounds
        % 0 if \pm\infty
        parameters.lambdam = beta-alpha;
        parameters.lambdap = beta+alpha;
        
        % F&L grid parameters
        parameters.FLc = delta;
        parameters.FLnu = 1; 
        
    case 3 % VG 
       %Schoutens page 57
       %C = 1.3574;
       %G = 5.8704;
       %M = 14.2699;
       %F-L
       C = 4;
       G = 12;
       M = 18;

        nu = 1/C;
        theta = (1/M-1/G)*C;
        s = sqrt(2*C/(M*G));

        % % Avramidis
        % nu = 0.4983;
        % theta = -0.28113;
        % s = 0.19071;

        % Rearrange parameters (time rescaling)
        parameters.nu = nu/dt;
        parameters.theta = theta*dt;
        parameters.s = s*sqrt(dt);

        % Compute vol(T), Schoutens page 60
        variance = s^2 + nu*theta^2;
        %parameters.vol = sqrt(variance*T);

        % analiticity bounds
        % 0 if +\infty
        parameters.lambdam = -M;
        parameters.lambdap = G;
        
    case 4 % Meixner, Schoutens page 62

        alpha = 0.3977;
        beta = -1.4940;
        delta = 0.3462;

        % Rearrange parameters (time rescaling)
        parameters.alpha = alpha;
        parameters.beta = beta;
        parameters.delta = delta*dt;

        % Compute vol(T), Schoutens page 60
        s = alpha*sqrt(delta)/(4*cos(beta/2));
        %parameters.vol = s*sqrt(T);

    case 5 % CGMY, Schoutens page 60

        %C = 0.0244; % C > 0
        %G = 0.0765; % G > 0
        %M = 7.5515; % M > 0
        %Y = 1.2945; % Y < 2
        %Y = 0.001; % Y < 2
        C = 4; 
        G = 50;
        M = 60;
        Y = 0.7;

        % Rearrange parameters (time rescaling)
        parameters.C = C*dt;
        parameters.G = G;
        parameters.M = M;
        parameters.Y = Y;

        % Compute vol(T), Schoutens page 60
        variance = C*(M^(Y-2)+G^(Y-2))*gamma(2-Y);
        %parameters.vol = sqrt(variance*T);
        
        % analiticity bounds
        % 0 if +\infty
        parameters.lambdam = -M;
        parameters.lambdap = G;
        
        % F&L grid parameters
        parameters.FLc = 2*C*abs(gamma(-Y)*cos(pi*Y/2));
        parameters.FLnu = Y;

    case 6 % Kou double exponential, Fusai and Roncoroni page 53

%         s = 0.120381;
%         lambda = 0.330966;
%         pigr = 0.20761; 
%         eta1 = 9.65997;
%         eta2 = 3.13868;
        s = 0.1;
        lambda =3;
        pigr = 0.3;
        eta1 = 40;
        eta2 = 12;
        
        % Rearrange parameters (time rescaling)
        parameters.s = s*sqrt(dt);
        parameters.lambda = lambda*dt;
        parameters.pigr = pigr;
        parameters.eta1 = eta1;
        parameters.eta2 = eta2;

        % Compute vol(T)
        variance = s^2 + lambda*(pigr/eta1^2+(1-pigr)/eta2^2);
        %parameters.vol = sqrt(variance*T);
        
        % analiticity bounds
        % 0 if +\infty
        parameters.lambdam = -eta1;
        parameters.lambdap = eta2;
        
        % F&L grid parameters
        parameters.FLc = s^2/2;
        parameters.FLnu = 2;

    case 7 % Merton jump-diffusion, Feng-Linetsky & Fusai and Roncoroni page 53

        s = 0.4; % 0.1; % 0.126349;
        alpha = 0.1; % -0.05; % -0.390078;
        lambda = 0.5; % 3; % 0.174814;
        delta = 0.15; % 0.086; % 0.338796;

        % Rearrange parameters (time rescaling)
        parameters.s = s*sqrt(dt);
        parameters.alpha = alpha;
        parameters.lambda = lambda*dt;
        parameters.delta = delta;

        % Compute vol(T)
        variance = s^2 + lambda*(delta^2+alpha^2);
        %parameters.vol = sqrt(variance*T);
        
        % analiticity bounds
        % 0 if +\infty
        parameters.lambdam = 0;
        parameters.lambdap = 0;
        
        % F&L grid parameters
        parameters.FLc = s^2/2;
        parameters.FLnu = 2;

    case 8 % Stable, Le Courtois and Quittard-Pinon DEF 31, 51-72 (2008)

        alpha = 2; % 1.4;
        beta = 0; % 0.236;
        gamm = 0.3/sqrt(2); % 0.15;
        m = 0;
        c = 1/sqrt(1+(beta*tan(alpha/2*pi))^2);

        % Rearrange parameters (time rescaling)
        parameters.alpha = alpha;
        parameters.beta = beta;
        parameters.gamm = gamm*dt;
        parameters.m = m;
        parameters.c = c;
        
        % Compute vol(T)
        if alpha == 2
            s = sqrt(2)*gamm;
        else
            s = Inf;
        end
        %parameters.vol = s*sqrt(T);

    otherwise

    e = 1;

end
