function [x,h,xi,H] = kernel(ngrid,xmin,xmax,parameters,alpha,disc,flag)
% disc = 0: no discount factor in the density
% disc = 1: discount factor in the density
% flag = 0: characteristic function for backward problem 
% flag = 1: characteristic function for forward problem

if nargin == 5
    disc = 1; % discount factor
    flag = 0; % backward characteristic function
elseif nargin == 6
    flag = 0; % backward characteristic function
end
    
N = ngrid/2;
dx = (xmax-xmin)/ngrid;
x = dx*(-N:N-1);
dxi = 2*pi/(xmax-xmin);
xi = dxi*(-N:N-1);
if nargin < 5 % shift parameter, esp. for Feng-Linetsky and convolution
    alpha = 0;
end
H = charfunction(xi+1i*alpha,parameters,flag); % characteristic function
if disc==1
    H = H*exp(-parameters.rf*parameters.dt); % discount
end
h = real(fftshift(fft(ifftshift(H))))/(xmax-xmin); % discounted kernel

% figure
% plot(xi,real(H),'r',xi,imag(H),'g')
% xlabel('xi')
% ylabel('\Psi(\xi,\Delta t)')
% legend('Re \Psi(\xi,\Delta t)','Im \Psi(\xi,\Delta t)')
% title('Characteristic function')
% 
% figure
% plot(x,h)
% xlabel('x')
% ylabel('f(x,\Delta t)')
% title('Probability density function')