% Jingbang code
close all
clear
tic

temp_seed = 4;

rng(temp_seed)
disp("seed is "+num2str(temp_seed))


Nx = 64;

steps = 1e15;
% rho0 = 0.8;
eps = 0.00011; % noise amplitude
everyData = 1; % how many time steps before checking for a rare event (checking costs time)
everyPlot = 5e5;  % hoe many time steps before plotting the shape

disp("eps is "+num2str(eps))

kappa = 0.008;
Crr = 6;
D0 = 1;
L = 1;
rho0_add = 0.01;
noiseon = 1;
% calculate linearly unstable rho0
ctilde = 1/(Crr-2*pi^2*kappa);
rho0_unstable = (1+sqrt(1-4*ctilde))/2;
rho0 = rho0_unstable + rho0_add;
format long
disp("rho0 is "+num2str(rho0))

x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;
dt = 1*500/((8*dk)^4);

Mob = rho0*(1-rho0)*D0*exp(-Crr*rho0);
Op = -Mob.*kappa/2.*k.^4;

N = @(rhoR) -(log(rhoR)-log(1-rhoR)-Crr.*rhoR); 

M1 = exp(Op*dt);
M2 = (exp(Op*dt)-1)./Op;
M2(1) = 0*dt;
M3 = sqrt(0.5*(exp(2*Op*dt)-1)./Op);
M3(1) = 0*sqrt(dt); 

E = @(rho) rho.*log(rho) + (1-rho).*log(1-rho) - 0.5*Crr*rho.^2;
Efunc = @(rho) sum(rho.*log(rho) + (1-rho).*log(1-rho) - 0.5*Crr*rho.^2)*dx;

% dV = @(rho,rhoR) ifft(-Op.*rho,"symmetric")+log(rhoR)-log(1-rhoR)-Crr.*rhoR;

time=0;
rhoRmin = 0.174726098877148;
Nevents = 25;  % number of events required

events = []; waiting_time = [];

rhoR = zeros(1,Nx)+rho0;
rho = fft(rhoR);
for step = 1:steps
    time = time + dt;
    xi = sqrt(2*Mob*eps/dx)*randn(size(rho));      % noise term 
    rhoR_ = rhoR;
    rho = M1.*rho + M2.*Mob.*(k.^2).*fft(N(rhoR)) + noiseon*M3.*(1i.*k).*fft(xi);
    rhoR = ifft(rho,"symmetric");
    if mod(step, everyPlot) == 0
%         clf; plot(x, rhoR); drawnow
%         clf; plot(k,abs(rho)); drawnow
%         pause
        disp("number of event is "+num2str(size(events,1)))
    end
    if mod(step,everyData)==0
        if any(rhoR<rhoRmin)                  % any event?
            [rhoRnow idx] = min(rhoR);        % if more than one event, pick min
            % interpolate linearly between current and previous step to
            % find the time of crossing
            rhoRthen = rhoR_(idx); t_int = (rhoRmin-rhoRthen)*dt/(rhoRnow-rhoRthen); 
            rhoR_int = rhoR_ + (rhoR-rhoR_)*t_int/dt;
            events = [events; circshift(rhoR_int,-idx+Nx/2+1)]; % store the event with lowest point at x=L/2
            rhoR = zeros(1,Nx)+rho0; rho = fft(rhoR); % reset rhoR for restart
            waiting_time = [waiting_time; time];
            time = 0;
        end
    end
    if size(events,1)>=Nevents
        break
    end
end

toc

% save workspace
filename = "rho0_"+num2str(rho0)+"_eps_"+num2str(eps)+"_"+num2str(temp_seed)+".mat";
save(filename)

exit
