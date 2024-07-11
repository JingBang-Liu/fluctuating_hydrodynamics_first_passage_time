clear
close all
Ncopies = 1e3;
maxsteps = 1e8;
dt = 5e-3;

% double well with varying mobility, and a zero eigen vector
alpha = 2;
beta = 1.5;
gamma = 1;
p1 = 1;
p2 = 1/4;
p = [p1;p2]; p = p./norm(p);
m = [-p(2);p(1)];
S_ = @(x, v) [sqrt(0.5.*(1+x(1,:).^2)).*p(1)^2.*v(1,:)+...
              sqrt(0.5.*(1+x(1,:).^2)).*p(1).*p(2).*v(2,:);...
              sqrt(0.5.*(1+x(1,:).^2)).*p(1).*p(2).*v(1,:)+...
              sqrt(0.5.*(1+x(1,:).^2)).*p(2)^2.*v(2,:)] ;
M_ = @(x, v) [0.5.*(1+x(1,:).^2).*p(1)^2.*v(1,:)+...
              0.5.*(1+x(1,:).^2).*p(1).*p(2).*v(2,:);...
              0.5.*(1+x(1,:).^2).*p(1).*p(2).*v(1,:)+...
              0.5.*(1+x(1,:).^2).*p(2)^2.*v(2,:)];
M = @(x) [0.5.*(1+x(1,:).^2).*p(1).^2,...
          0.5.*(1+x(1,:).^2).*p(1).*p(2);...
          0.5.*(1+x(1,:).^2).*p(1).*p(2),...
          0.5.*(1+x(1,:).^2).*p(2).^2];
divM = @(x) [x(1,:).*p(1).^2;...
             x(1,:).*p(1).*p(2)];

U = @(x) 0.25*(1-x(1,:).^2).^2 + x(2,:).^2./2.*(x(1,:).^2+0.25);
dU = @(x) [-x(1,:)+x(1,:).^3 + x(2,:).^2.*x(1,:);...
           x(2,:).*(x(1,:).^2+0.25)];
b = @(x) -M_(x, dU(x));
HessU = @(x) [3*x(1,:).^2-1+x(2,:).^2, 2.*x(1,:).*x(2,:);...
              2.*x(1,:).*x(2,:), x(1,:).^2+0.25];

% local minima, calculated by mathematica
a = [-0.9730518486301645;...
     -0.12144300123917997];

% saddle, calculated by mathematica
s = [0.008577885905028669;...
     0.1310534692073233];


temp = sort(eig(HessU(s)));
lambda_minus = temp(1);
temp2 = sort(eig(HessU(s)*M(s)));
mu_minus = temp2(1);
[tempV,tempD] = eig(HessU(s)*M(s));

prefactor_naive = pi/abs(lambda_minus)*sqrt(-det(HessU(s))/det(HessU(a)));
prefactor_mobility = pi/abs(mu_minus)*sqrt(-det(HessU(s))/det(HessU(a)));
prefactor_mobility_cons = pi/abs(mu_minus)*sqrt(-det(HessU(s))*(m'*inv(HessU(s))...
    *m)/det(HessU(a))/(m'*inv(HessU(a))*m));
deltaU = U(s)-U(a);

mfpt_naive = @(eps) prefactor_naive*exp(eps.^-1*deltaU);
mfpt_mobility = @(eps) prefactor_mobility*exp(eps.^-1*deltaU);
mfpt_mobility_cons = @(eps) prefactor_mobility_cons*exp(eps.^-1*deltaU);

epsis = 10.^linspace(0.5,log(1e-1)/log(10),15);
mfpts = 0*epsis;
figure(1)
for idx = 1:size(epsis,2)
  eps = epsis(idx);
  b = @(x) -M_(x, dU(x)) + eps*divM(x);
  check_prefactor_sampling_varying_0
  mfpts(idx) = mean(times)

  clf;
  semilogy(1./epsis, mfpt_naive(epsis),'b--'); hold on
  semilogy(1./epsis, mfpt_mobility(epsis),'g--'); hold on
  semilogy(1./epsis, mfpt_mobility_cons(epsis),'k--'); hold on
  semilogy(1./epsis, mfpts, 'rx')
  title('MFPT double well, varying mobility, conserved quantity', 'interpreter', 'latex')
  legend({'Original Eyring-Kramers',...
		  'With mobility correction',...
          'With conserved quantity','sampling'}, 'location',...
		 'best', 'interpreter', 'latex')
  xlabel('$1/\varepsilon$', 'interpreter', 'latex')
  grid
  drawnow
end


