function figures_ft(S,x,xi,f,Psi,g,G)

figure
plot(xi,real(Psi),'r',xi,imag(Psi),'g')
xlabel('xi')
ylabel('\Psi(\xi,T)')
legend('Re \Psi(\xi,T)','Im \Psi(\xi,T)')
title('Characteristic function at time T')

figure
plot(x,real(f),'r',x,imag(f),'g')
xlabel('x')
ylabel('f(x,T)')
title('Probability density function at time T')

% Comparison of the numerical and analytical Fourier transforms of the payoff

% Normal space
gn = fftshift(fft(ifftshift(G)))/((x(2)-x(1))*length(x));

figure
plot(x,g,'r',x,real(gn),'g')
xlabel('x')
ylabel('Re g')
legend('analytical','numerical')
title('Payoff function')

figure
plot(x,zeros(size(x)),'ro',x,imag(gn),'gs')
xlabel('x')
ylabel('Im g')
legend('analytical','numerical')
title('Payoff function')

figure
plot(x,g./real(gn),'gs')
xlabel('x')
ylabel('g_a/Re g_n')
xlim([0 u])
title('Payoff function')

figure
plot(x,g,'r',S,real(gn),':r',x,imag(gn),':g')
xlabel('x')
ylabel('g')
legend('g_a','Re g_n','Im g_n')
title('Payoff function')

% Reciprocal space
Gn = fftshift(ifft(ifftshift(g)))*(x(2)-x(1))*length(x);

figure
plot(xi,real(G),'ro',xi,real(Gn),'gs')
%xlim([-50 50])
xlabel('xi')
ylabel('Re G')
legend('analytical','numerical')
title('Fourier transform of the payoff function')

figure
plot(xi,imag(G),'ro',xi,imag(Gn),'gs')
%xlim([-50 50])
xlabel('xi')
ylabel('Im G')
legend('analytical','numerical')
title('Fourier transform of the payoff function')

figure
plot(xi,real(G)./real(Gn),'ro',xi,imag(G)./imag(Gn),'gs')
%xlim([-50 50])
xlabel('xi')
legend('Re(G_a)/Re(G_n)','Im(G_a)/Im(G_n)')
title('Fourier transform of the payoff function')

figure
plot(xi,real(G),'r',xi,real(Gn),':r',xi,imag(G),'g',xi,imag(Gn),':g')
%xlim([-50 50])
xlabel('\xi')
ylabel('G')
legend('Re G_a','Re G_n)','Im G_a','Im G_n')
title('Fourier transform of the payoff function')