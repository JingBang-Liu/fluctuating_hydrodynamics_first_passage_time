close all
clear


%% calculate prefactors
Crr = 6;
kappa = 0.008;
D0 = 1;
L = 1;
Nx = 128;
rho0_add = 0.01;
ctilde = 1/(Crr-2*pi^2*kappa);
rho0_unstable = (1+sqrt(1-4*ctilde))/2;
rho0 = rho0_unstable + rho0_add;
Mob = rho0*(1-rho0)*D0*exp(-Crr*rho0);


x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;

% GAD_SEG_JB_for_prefactor
load("s_128.mat")

a = 0*x+mean(s);

E = @(rhoR) mean(rhoR.*log(rhoR)+(1-rhoR).*log(1-rhoR)-0.5*Crr.*rhoR.^2 ...
    -0.25*kappa.*rhoR.*ifft(-k.^2.*fft(rhoR),"symmetric"))*L;

Us = E(s);
Ua = E(a);
deltaU = Us - Ua;
disp("deltaU "+num2str(deltaU))

M = -Mob*ifft(-k'.^2.*fft(eye(Nx)), 'symmetric');
Ms = -ifft(1i*k'.*fft(diag(s.*(1-s)*D0.*...
    exp(-Crr.*s-0.5*kappa.*ifft(-k.^2.*fft(s))))*ifft(1j*k'.*fft(eye(Nx)),"symmetric")));
Hs = -(-ifft(0.5*kappa.*(-k'.^2).*fft(eye(Nx)) + fft(diag(s.^(-1)+(1-s).^(-1)-Crr)), 'symmetric'));
Ha = -(-ifft(0.5*kappa.*(-k'.^2).*fft(eye(Nx)) + fft(diag(a.^(-1)+(1-a).^(-1)-Crr)), 'symmetric'));
[VHMs,DHMs] = eig(Hs*M);
[VHMss,DHMss] = eig(Hs*Ms);
[VMHss,DMHss] = eig(Ms*Hs);
[VMHs,DMHs] = eig(M*Hs);
[VHs,DHs] = eig(Hs);
[VHa,DHa] = eig(Ha);
[VM,DM] = eig(M);
[VMs,DMs] = eig(Ms);
spec_HMs = zeros(Nx,1);
spec_HMss = zeros(Nx,1);
spec_MHss = zeros(Nx,1);
spec_MHs = zeros(Nx,1);
spec_Hs = zeros(Nx,1);
spec_Ha = zeros(Nx,1);
spec_M = zeros(Nx,1);
spec_Ms = zeros(Nx,1);
for n=1:Nx
    spec_HMs(n) = real(DHMs(n,n));
    spec_HMss(n) = real(DHMss(n,n));
    spec_MHss(n) = real(DMHss(n,n));
    spec_MHs(n) = real(DMHs(n,n));
    spec_Hs(n) = real(DHs(n,n));
    spec_Ha(n) = real(DHa(n,n));
    spec_M(n) = real(DM(n,n));
    spec_Ms(n) = real(DMs(n,n));
end
[spec_HMs,IHMs] = sort(spec_HMs,'descend');
[spec_HMss,IHMss] = sort(spec_HMss,'descend');
[spec_MHss,IMHss] = sort(spec_MHss,'descend');
[spec_MHs,IMHs] = sort(spec_MHs,'descend');
[spec_Hs,IHs] = sort(spec_Hs,'descend');
[spec_Ha,IHa] = sort(spec_Ha,'descend');
[spec_M,IM] = sort(spec_M);
[spec_Ms,IMs] = sort(spec_Ms);
VHMs(:,:) = VHMs(:,IHMs);
VHMss(:,:) = VHMss(:,IHMss);
VMHss(:,:) = VMHss(:,IMHss);
VMHs(:,:) = VMHs(:,IMHs);
VHs(:,:) = VHs(:,IHs);
VHa(:,:) = VHa(:,IHa);
VM(:,:) = VM(:,IM);
VMs(:,:) = VMs(:,IMs);

%%% prefactor
mu = abs(spec_HMs(1)+spec_HMs(2));
mu_s = abs(spec_HMss(1)+spec_HMss(2));
lambda = abs(spec_Hs(1));
nhat = VM(:,1);
nhat_s = VMs(:,1);
nHsn = abs(nhat'*inv(Hs)*nhat);
nHan = abs(nhat'*inv(Ha)*nhat);
nHsn_s = abs(nhat_s'*inv(Hs)*nhat_s);
nHan_s = abs(nhat_s'*inv(Ha)*nhat_s);
prefactor_mob_var_mass = 2*pi/mu_s*sqrt(abs(prod(spec_Hs./spec_Ha)))*sqrt(nHsn_s/nHan_s); % conserved quantity and variable mobility
prefactor_mass = 2*pi/mu*sqrt(abs(prod(spec_Hs./spec_Ha)))*sqrt(nHsn/nHan); % with conserved quantity
prefactor_mob = 2*pi/mu*sqrt(abs(prod(spec_Hs./spec_Ha))); % with mobility
prefactor_original = 2*pi/lambda*sqrt(abs(prod(spec_Hs./spec_Ha))); % original
disp("prefactor original "+num2str(prefactor_original))
disp("prefactor mobility "+num2str(prefactor_mob))
disp("prefactor mass     "+num2str(prefactor_mass))
disp("prefactor mob var  "+num2str(prefactor_mob_var_mass))


 %% calculate time
% epsis = 0.34*10.^linspace(0,log(1e-1)/log(10),15);
% epsis = [0.1,0.01,0.001,0.0004,0.0003,0.0002,0.0001,0.00009,0.00008,0.00007,0.00006,0.00005,0.00004,0.00003,0.00002];
epsis = [0.1,0.01,0.001,0.0004,0.0003,0.0002,0.00016,0.00014,...
    0.00012,0.00011,0.0001,0.00009,0.00008,0.00007,0.00006,0.00005,0.00004,0.00003,0.000025];


tau_Kramar_mob_var_mass= prefactor_mob_var_mass.*exp(deltaU./epsis);
tau_Kramar_mob_mass= prefactor_mass.*exp(deltaU./epsis);
tau_Kramar_mob = prefactor_mob.*exp(deltaU./epsis);
tau_Kramar_original = prefactor_original.*exp(deltaU./epsis);

epsis_MC = [0.1,0.01,0.001,0.0004,0.0003,0.0002,0.00016,0.00014,0.00012,0.00011,0.0001,0.00009,0.00008,0.00007,0.00006,0.00005,0.00004,0.00003];
% tau_MC = [8.1065e-04,0.1486,85.1218,1.2577935e+02,1.5380301e+02,...
%     1.934732e+02,313.3028,3.130642e+02,3.73562e+02,4.513272e+02,...
%     5.317408e+02,641.8335,8.0184613e+02,1.5947822e+03,4.5247e+03];
tau_MC = [8.1065e-04,0.1486,85.1218,1.2577935e+02,1.5380301e+02,...
    1.934732e+02,2.227634e+02,2.31483e+02,2.5799201e+02,3.0715743e+02,...
    313.3028,3.130642e+02,3.73562e+02,4.513272e+02,...
    5.317408e+02,641.8335,8.0184613e+02,1.5947822e+03];

%% plot
semilogy(1./epsis, tau_Kramar_original,"Color", [0.3010 0.7450 0.9330],"LineStyle","--"); hold on
semilogy(1./epsis, tau_Kramar_mob,"Color",[0 0.4470 0.7410],"LineStyle","--"); hold on
semilogy(1./epsis, tau_Kramar_mob_mass,'r-'); hold on
semilogy(1./epsis, tau_Kramar_mob_var_mass,'g-'); hold on
semilogy(1./epsis_MC, tau_MC,"k."); hold on
% semilogy(1./epsis_1, tau_MC_1,"r*"); hold on
title('Eyring-Kramars, rho0 = rho0Stable+0.01', 'interpreter', 'latex','fontsize',16)
legend({'original Eyring-Kramars',...
		'with mobility correction', ...
		'with conserved quantity',...
        'with non-constant mobility',...
        'Monte-Carlo simulation const-mob'}, 'location',...
	   'best', 'interpreter', 'latex', 'fontsize', 16)
ylabel('mean first passage time $w_B(x_)$', 'interpreter', 'latex',"FontSize",15)
xlabel('noise strength $1/\varepsilon$', 'interpreter', 'latex',"FontSize",15)
grid

% now, save the data for python plot script
save_data_theory = zeros(size(epsis,2),4);
save_data_theory(:,1) = epsis;
save_data_theory(:,2) = tau_Kramar_original;
save_data_theory(:,3) = tau_Kramar_mob;
save_data_theory(:,4) = tau_Kramar_mob_mass;
writematrix(save_data_theory,"SEG_epsis_tau_theory.txt",'Delimiter','space')
save_data_MC = zeros(size(epsis_MC,2),2);
save_data_MC(:,1) = epsis_MC;
save_data_MC(:,2) = tau_MC;
writematrix(save_data_MC,"SEG_epsis_tau_MC.txt",'Delimiter','space')

