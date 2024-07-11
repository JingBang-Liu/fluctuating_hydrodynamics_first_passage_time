% integrates d_t h = h0^3 d_xxxx h + sqrt(2 eps h0^3) d_x eta
% for spatio-temporal white eta
clear
tic
h0 = 1.01;
steps = 1e15;

temp_seed = 1;

rng(temp_seed)
disp("seed is "+num2str(temp_seed))

startup = 0;  %0 from flat film, 1 from thermal equilibrium  

z_hydro = 0; %Old shit.

L = 1;       % domain size
eps = 0.0034;  % noise strength
noiseon=1;   %Turn noise on or off
ell=0;

disp("eps is "+num2str(eps))

std_ana = sqrt(eps/12);

Nx = 128;    % number of nodes

Nevents = 200  % number of events required
z = z_hydro;   %0.874*h0; % z     %h0=1.05->0.777;  1.26->0.556  % minimum height which corresponds to a rare event - z=0 for breakup, but z>0 to find saddle points


for i=1:1

  Mob = (h0-z_hydro)^3+ell*h0^2;

  Ncopies = 1; % this gives the ability to run multiple realisations simultaneously, which could be faster (but doesn't seem to be!).

  everyData = 1; % how many time steps before checking for a rare event (checking costs time)
  everyPlot = 1e10;  % hoe many time steps before plotting the shape

  % Spatial and fourier grids
  x = linspace(0,L-L/Nx,Nx);
  dx = x(2)-x(1);
  dk = 2*pi/L;
  k = [0:Nx/2,-Nx/2+1:-1]*dk;

  % Relaxation of modes (from FT of d_xxxx h )
  c = -(Mob)*k.^4;

  % This is the tricky part.  How to choose dt?  So that all modes are
  % captured (painful), or only some?  Here you see one choice..??????

  if(h0<=0.6)
    dt = 1*0.001/((8*dk)^4); %Resolve the first four modes
  elseif(h0<=0.8)
    dt = 1*0.01/((8*dk)^4); %Resolve the first four modes
  elseif(h0<=1.0)
    dt = 1*0.1/((8*dk)^4); %Resolve the first four modes
  else
    dt = 1*1/((8*dk)^4); %Resolve the first four modes
  end

  %dt = 1*0.1/((8*dk)^4) %Resolve the first four modes

  % setting up first order exponential time differences
  M1 = exp(c*dt);
  M1 = repmat(M1,[Ncopies,1]);
  M2 = (exp(c*dt)-1)./c;
  M2(1) = 0*dt;

  % this makes copies so you can run multiple realisations at once
  M2 = repmat(M2,[Ncopies,1]);

  M3 = sqrt(0.5*(exp(2*c*dt)-1)./c);

  M3(1) = 0*sqrt(dt); % do not force mean
  M3 = repmat(M3,[Ncopies,1]);

  % a function calculating the disjoining pressure term in real space
  N = @(h) -Mob*(4*pi^2/3)*(ifft(h,[],2, 'symmetric')).^-3;

  events = []; waiting_time = [];

  % note 'h' is in Fourier space and hR is in real space
  if(startup==0)
    hR = zeros(Ncopies,Nx) + h0;
    h = fft(hR,[],2);                                   
  else
    h=FTh0(h0,eps,x,Nx);
    hR = ifft(h,[],2, 'symmetric');
  end

  % here I just start from a flat interface, but I can easily enough start
  % from realisations drawn from thermal equilibrium 

  time(Ncopies,1)=0;

  for step = 1:steps
    time(:,1) = time(:,1)+dt;  hR_=hR;
    
    eta = sqrt(2*Mob*eps/dx)*randn(size(h));      % noise term 
    h = M1.*h + M2.*(k.^2).*fft(N(h),[],2) + noiseon*M3.*(1i.*k).*fft(eta,[],2);   % time step in Fourier space
    
    % checking for rare events
    if mod(step,everyData)==0
      if mod(step, 1e7)==0
        disp("this step "+num2str(step)+" with "+num2str(size(events,1))+" ruptures")
      end           
      hR = ifft(h,[],2, 'symmetric'); 

      if any(any(hR<z))                    % any event?
        for c = 1:Ncopies                  % run through realisations
          idxs = find(hR(c,:)<z);          % find position if there is an event
          if(size(idxs,2)==0)               % skip this c is there is no event  
            continue
          else
            [hnow idx] = min(hR(c,:));           % if there is more than one point satisfying the rare event critetion choose the minimum
          end
            
          %interpolate linearly between current and previous time step to
          %find time at which curves crosses the rare event boundary
          %this is potentially dodgy and could be corrected with smaller time
          %steps / Brownian bridges.
          hnow=hR(c,idx); hthen=hR_(c,idx); t_int = (z-hthen)*dt/(hnow-hthen);
          hR_int(c,:) = hR_ + (hR-hR_)*t_int/dt;

          events = [events; circshift(hR_int(c,:), -idx+Nx/2+1)];  % store the event with lowest point at x=L/2
          if(startup==0)
            h(c,:) = fft(zeros(1,Nx) + h0,[],2);                                       % reset the realisation with an event
          else
            h(c,:) = FTh0(h0,eps,x,Nx);
          end
          waiting_time = [waiting_time; time(c,1)];                                  % store waiting time
          time(c,1) = 0;                                                             % reset time
        end
      end
      if size(events,1)>=Nevents
        break
      end
    end
  end

  waiting_time_is(i)=mean(waiting_time)
end
toc



% save workspace
filename = "h01.01_"+num2str(eps)+"_"+num2str(temp_seed)+".mat";
save(filename)

exit
