function F = charfunction(xi,parameters,flag)
% flag=0 --> funzione caratteristica per problema backward 
% flag=1 --> funzione caratteristica per problema forward

if nargin==2
    flag=0;
end

meancorrection = (parameters.rf-parameters.q)*parameters.dt-log(charfunction0(-1i,parameters));
F = exp(1i*meancorrection*xi).*charfunction0(xi,parameters);
if flag==0
    F=conj(F);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = charfunction0(xi,parameters)

switch parameters.distr

    case 1 % Normal
        m = parameters.m;
        s = parameters.s;
        F = exp(1i*xi*m-0.5*(s*xi).^2);

    case 2 % Normal inverse Gaussian (NIG) by Barndorff-Nielsen (1998), Schoutens page 59
        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta;
        F = exp(-delta*(sqrt(alpha^2-(beta+1i*xi).^2)-sqrt(alpha^2-beta^2)));

    case 3 % Variance Gamma (VG) by Madan and Seneta (1990), Schoutens page 57
        theta = parameters.theta;
        s = parameters.s;
        nu = parameters.nu;
        F = (1-1i*xi*theta*nu+0.5*nu*(s*xi).^2).^(-1/nu);

    case 4 % Meixner, Schoutens page 62
        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta;
        F = (cos(beta/2)./cosh((alpha*xi-1i*beta)/2)).^(2*delta);

    case 5 % Carr, Geman, Madan and Yor (CGMY) model (2002), Schoutens page 60
        C = parameters.C;
        G = parameters.G;
        M = parameters.M;
        Y = parameters.Y;
        F = exp(C*gamma(-Y)*((M-1i*xi).^Y-M^Y+(G+1i*xi).^Y-G^Y));

    case 6 % Kou double exponential jump-diffusion (2002), Fusai and Roncoroni page 53
        s = parameters.s;
        lambda = parameters.lambda;
        pigr = parameters.pigr;
        eta1 = parameters.eta1;
        eta2 = parameters.eta2;
        F = exp(-0.5*(s*xi).^2+lambda*((1-pigr)*eta2./(eta2+1i*xi)+pigr*eta1./(eta1-1i*xi)-1));

    case 7 % Merton jump-diffusion (1976), Fusai and Roncoroni page 53
        s = parameters.s;
        alpha = parameters.alpha;
        lambda = parameters.lambda;
        delta = parameters.delta;
        F = exp(-0.5*(s*xi).^2+lambda*(exp(1i*xi*alpha-0.5*(delta*xi).^2)-1));

    case 8 % Levy alpha-stable
        alpha = parameters.alpha;
        beta = parameters.beta;
        gamm = parameters.gamm;
        m = parameters.m;
        c = parameters.c;
        F = exp(1i*xi*m-c*abs(gamm*xi).^alpha.*(1-1i*beta*sign(xi)*tan(alpha/2*pi)));

    otherwise
    e = 1;

end
