% second derivative, not fourth derivative, so only v, no w

dt = 0.001;
steps = 2e6;
everyPlot = 1e4;

%%% some parameters
h0 = 1.01;
ell = 0.0;
L = 1;
Nx = 128;
magnitude = 0.1;
Mob = 1;

x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;

Op = -k.^2;

% setting up first order exponential time differences
M1 = exp(Op*dt);
M2 = (exp(Op*dt)-1)./Op;
M2(1) = 0*dt;

hR = h0 + magnitude*cos(2*pi*x/L);
h = fft(hR);
vR = cos(2*pi*x/L);
vR = vR./sqrt(sum(vR.*vR));
v = fft(vR);
v(1) = 0;
figure(1)
for step = 1:steps
  h_ = h;
  hR = ifft(h, 'symmetric');
  vR = ifft(v, 'symmetric');
  vR = vR./sqrt(sum(vR.*vR));
  v = fft(vR);
  dV = ifft(k.^2.*h, 'symmetric') + 4*pi^2/3./(hR.^3);
  ddVv = ifft(k.^2.*v, 'symmetric')-4*pi^2./(hR.^4).*vR;
  Cx = sum(dV.*vR)/sum(vR.*vR);
  Calpha = sum(ddVv.*vR);
  h = M1.*h + M2.*fft(-4*pi^2/3./hR.^3 + 2*Cx*vR);
  v = M1.*v + M2.*fft(4*pi^2./hR.^4.*vR + Calpha*vR);
  if mod(step, everyPlot) == 0
	clf; plot(x, hR); drawnow
    disp(max(Mob*abs(ifft(k.^4.*fft(hR)+4*pi^2/3*k.^2.*fft(1./hR.^3),'symmetric'))))
  end
  hR = ifft(h,'symmetric');
  % make sure change is small and dh/dt is small
  if (max(abs(ifft(h)-ifft(h_)))<1e-15)&& ...
      (max(Mob*abs(ifft(k.^4.*fft(hR)+4*pi^2/3*k.^2.*fft(1./hR.^3),'symmetric')))<1e-4)
    disp("break")
    break
  end
end

clf; 
plot(x, hR)
ylim([0,max(hR)+0.2])

hSaddle = hR;
s = hSaddle;
save("s_temp_1.01.mat","s");
% close all
