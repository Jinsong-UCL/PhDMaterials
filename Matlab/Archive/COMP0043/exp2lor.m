%% Check numerically the Fourier pair Laplace <-> Lorentzian

% Grids in real and Fourier space
% They are linked by the Nyquist relation Dx*Dxi = 2*pi/N
N = 2048; % grid size
Dx = 0.01; % grid step in real space
Lx = N*Dx; % upper truncation limit in real space
Dxi = 2*pi/Lx; % grid step in Fourier space
Lxi = N*Dxi; % upper truncation limit in Fourier space
x = Dx*(-N/2:N/2-1); % grid in real space
xi = Dxi*(-N/2:N/2-1); % grid in Fourier space

% Analytical expressions
a = 1;
fa = a/2*exp(-a*abs(x)); % Laplace (or double exponential)
Fa = a^2./(a^2+xi.^2); % Lorentzian (or Cauchy)

Fn  = fftshift(ifft(ifftshift(fa)))*Lx;
fn  = fftshift( fft(ifftshift(Fa)))/Lx;
Fn1 = fftshift( fft(ifftshift(fa)))*Dxi;
fn1 = fftshift(ifft(ifftshift(Fa)))/Dxi;

close all
figure(1), clf, hold on
plot(x,real(fn),'r')
plot(x,imag(fn),'g')
plot(x,fa,'k:')
axis([-10 10 0 1])
xlabel('x')
ylabel('f')
legend('Re(fn)','Im(fn)','fa')

figure(2), clf, hold on
plot(xi,real(Fn),'b')
plot(xi,imag(Fn),'m')
plot(xi,Fa,'c:')
axis([-20 20 0 1])
xlabel('\xi')
ylabel('F')
legend('Re(Fn)','Im(Fn)','Fa')

figure(3), clf, hold on
plot(x,real(fn1),'r')
plot(x,imag(fn1),'g')
plot(x,fa,'k:')
axis([-10 10 0 1])
xlabel('x')
ylabel('f')
legend('Re(fn1)','Im(fn1)','fa')

figure(4), clf, hold on
plot(xi,real(Fn1),'b')
plot(xi,imag(Fn1),'m')
plot(xi,Fa,'c:')
axis([-20 20 0 1])
xlabel('\xi')
ylabel('F')
legend('Re(Fn1)','Im(Fn1)','Fa')