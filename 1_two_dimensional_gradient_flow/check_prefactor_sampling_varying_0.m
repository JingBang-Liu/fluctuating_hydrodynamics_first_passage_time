% sampling only
x = zeros(2,Ncopies)+a;
times = zeros(Ncopies,1);

for step = 1:maxsteps
  %x = x + dt*b(x) + sqrt(2*dt*eps)*S(x)*randn(size(x));
  x = rk4_tv_step_varying( x, 0, dt, 2*eps, b, S_);
  if all(times>0)
	break
  end
  times(logical((times==0).*(x(1,:)'>s(1))))=step*dt;
end