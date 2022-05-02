N = 100;
s00 = S0*ones(N,1,'gpuArray');

finalStockPrices = arrayfun( @MC_GPU,nsteps,d,v_0,dt,rho_v,r_0,rho_r,a_i,a_j,b_i,b_j,param_alpha,s00,K)
