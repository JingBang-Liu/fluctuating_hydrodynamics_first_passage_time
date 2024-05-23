% GAD for MMM with constant mobility, using simplified GAD
% no need to include mobility, according to Tobias
close all
clear

steps = 2e8;
everyPlot = 1e4;

%%% some parameters
Crr = 6;
kappa = -0.008;
L = 1;
Nx = 128;
rho0_add = 0.01;
ctilde = 1/(Crr-2*pi^2*kappa);
rho0_unstable = (1+sqrt(1-4*ctilde))/2;
rho0 = rho0_unstable + rho0_add;
magnitude = rho0/10;

x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;

dt = 0.0001;

Op = kappa/2.*k.^2;

E = @(rhoR) mean(rhoR.*log(rhoR)+(1-rhoR).*log(1-rhoR)-0.5*Crr.*rhoR.^2 ...
    +0.25*kappa.*rhoR.*ifft(-k.^2.*fft(rhoR),"symmetric"))*L;

% setting up first order exponential time differences
M1 = exp(Op.*dt);
M2 = (exp(Op.*dt)-1)./Op;
M2(1) = 0*dt;

sigma = 0.1;
mu = L/2;
rhoR = rho0 + magnitude*cos(2*pi*x/L);
rho = fft(rhoR);
vR = cos(2*pi*x/L);
vR = vR./sqrt(sum(vR.*vR));
v = fft(vR);
figure
for step = 1:steps
  rho_ = rho;
  rhoR = ifft(rho, 'symmetric');
  vR = ifft(v, 'symmetric');
  v = fft(vR);
  dV = ifft(-Op.*rho,"symmetric")+log(rhoR)-log(1-rhoR)-Crr.*rhoR; % CH
  ddVv = ifft(-Op.*v,"symmetric")+(1./rhoR+1./(1-rhoR)-Crr).*vR;
  Cx = sum(dV.*vR)/sum(vR.*vR);
  Calpha = sum(ddVv.*vR)/sum(vR.*vR);
  Gx = -fft(log(rhoR)-log(1-rhoR)-Crr.*rhoR);
  Gv = -fft((1./rhoR+1./(1-rhoR)-Crr).*vR);
  rho = M1.*rho + M2.*(Gx + 2*Cx*v);
  v = M1.*v + M2.*(Gv + Calpha*v);

  if mod(step, everyPlot) == 0
	%clf; plot(x, phiR, x, vR); drawnow
	clf; plot(x, rhoR); drawnow
%     disp(mean(rhoR))
%     disp(max(abs(ifft(rho,"symmetric")-ifft(rho_,"symmetric"))))
    disp(E(rhoR))
%       disp(max(abs(dV)))
  end

  rhoR = ifft(rho,'symmetric');
  % make sure change is small and dh/dt is small
  if (max(abs(ifft(rho,"symmetric")-ifft(rho_,"symmetric")))<1e-11)&& ...
      (max(abs(ifft(-k.^2.*fft(dV),"symmetric")))<1e-9)
    disp("break")
    break
  end
end

% clf; 
% plot(x, rhoR)
% ylim([0,max(rhoR)+0.2])
% 
rhoSaddle = rhoR;
s = rhoSaddle;
save("s_128.mat","s")
close all
