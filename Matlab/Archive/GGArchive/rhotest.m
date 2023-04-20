N = 2;
%M = N;
M = 2*10000;

%rho = 2*rand(N)-1 % run a few times until positive semidefinite
rho = [-0.5417 0.1899;-0.1170 -0.4834]
NW = randn(M,N);
NB = randn(M,N);
%corr(NW,NB) % \approx zeros(N) if M >> N
%corr(NW,NW) % \approx eye(N)   if M >> N
%corr(NB,NB) % \approx eye(N)   if M >> N

NZ = NW*rho+NB*chol(eye(N)-rho.'*rho,'lower');
corr(NW,NZ) % \approx rho      if M >> N

rho_1 = zeros(N,N,10000);
for i = 1:10000
    N_W = NW(2*i-1:2*i,:);
    N_Z = NZ(2*i-1:2*i,:);
    rho_1(:,:,i) = N_W*N_Z.';
    %N_W*N_Z;
end
mean(rho_1,3)

