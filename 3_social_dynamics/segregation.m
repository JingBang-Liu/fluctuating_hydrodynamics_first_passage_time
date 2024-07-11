% Tobias code
close all
clear

Nx = 64;

dt = 0.001;
steps = 2e6;
everyPlot = 1e2;

kappa = 0.008;
Crr = 6;
L = 1;
rho0_add = 0.01; % how much added to rho0_unstable

x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;

Mob = k~=0;
Op = -Mob.*kappa/2.*k.^2;

N = @(rho) -log(rho) + log(1-rho)+ Crr*rho; 

M1 = exp(Op*dt);
M2 = (exp(Op*dt)-1)./Op;
M2(1) = 0*dt;

E = @(rho) rho.*log(rho) + (1-rho).*log(1-rho) - 0.5*Crr*rho.^2 - ...
    0.25*kappa.*rho.*ifft(-k.^2.*fft(rho),"symmetric");
Efunc = @(rho) sum(rho.*log(rho) + (1-rho).*log(1-rho) - 0.5*Crr*rho.^2 - ...
    0.25*kappa.*rho.*ifft(-k.^2.*fft(rho),"symmetric"))*dx;

dV = @(rhoR) -0.5*kappa.*ifft(-k.^2.*fft(rhoR),"symmetric")+log(rhoR)-log(1-rhoR)-Crr.*rhoR;

ctilde = 1/(Crr-2*pi^2*kappa);
rho0_unstable = (1+sqrt(1-4*ctilde))/2;
rho0 = rho0_unstable + rho0_add;
noise = 2e-1*randn(1,Nx);
noise = noise-mean(noise);
rhoR = zeros(1,Nx)+rho0+noise;
alpha = 0.6; delta = 0.16;
rhoR(1:floor(alpha*Nx)) = rho0 + delta;
rhoR(floor(alpha*Nx)+1:end) = rho0 - floor(alpha*Nx)/(Nx-floor(alpha*Nx))*delta;
mean(rhoR)

rho1R = 0*rhoR + rho0 + noise;
rho2R = rhoR;
rho1 = fft(rho1R);
rho2 = fft(rho2R);
for step = 1:steps
  rho1R_ = rho1R;
  rho2R_ = rho2R;
  rho1 = M1.*rho1 + M2.*(fft(N(rho1R)));
  rho2 = M1.*rho2 + M2.*(fft(N(rho2R)));
  rho1R = ifft(rho1, 'symmetric');
  rho2R = ifft(rho2, 'symmetric');
  if mod(step, everyPlot) == 0
	clf; plot(x, rho1R,x, rho2R); drawnow
  end
  if (max(abs(rho1R-rho1R_))<1e-13)&& ...
      (max(abs(rho2R-rho2R_))<1e-13)
    disp("break")
    break
  end
end

[mm,idx] = min(rho2R);
disp("minimum of rho is "+num2str(mm))
rho2R = circshift(rho2R,Nx/2-idx);

