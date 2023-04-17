N = 2;
rho = [-0.5417 0.1899;-0.1170 -0.4834]
%rho =  [-0.5417 0;0 -0.4834]
rho_1 = zeros(N,N,10000);
for i = 1:10000
    N_W = randn(N);
    N_B = randn(N);

    %N_Z = N_W*rho + N_B*chol(eye(N)-rho*rho.',"lower"); %Eq. (2.6) in Gnoatto's
    %rho_1(:,:,i) = N_W*N_Z.';
    %rho_1(:,:,i) = corr(NW,NZ);

    N_Z = rho*N_W.' + sqrt(eye(d)-rho.'*rho)*N_B.';% Eq. (167) in our paper
    rho_1(:,:,i) = corr(NW,NZ);
    %N_W*N_Z;
end
mean(rho_1,3)