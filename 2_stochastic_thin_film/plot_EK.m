close all
clear


%% calculate prefactors
Nx = 128;
h0 = 1.01;
ell = 0.0;
Mob = h0^3 + ell*h0^2;
L = 1;

x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;

load("s_temp_1.01.mat"); % load saddle shape
    
a = 0*x+mean(s);

Us = sum(1/2.*(ifft(1i*k.*fft(s),'symmetric')).^2-2/3*pi^2./(s.^2))*dx;
Ua = sum(1/2.*(ifft(1i*k.*fft(a),'symmetric')).^2-2/3*pi^2./(a.^2))*dx;
deltaU = (Us - Ua);
disp("deltaU "+num2str(deltaU))

M = -Mob*ifft(-k'.^2.*fft(eye(Nx)), 'symmetric');
Ms = -ifft( 1j*k'.*fft(diag(s.^3+ell.*s.^2)*ifft(1j*k'.*fft(eye(Nx)),"symmetric")));
Hs = ifft(-k'.^2.*fft(eye(Nx)) + 4*pi^2*fft(diag(s'.^(-4))), 'symmetric');
Ha = ifft(-k'.^2.*fft(eye(Nx)) + 4*pi^2*fft(diag(a'.^(-4))), 'symmetric');

[VHMs,DHMs] = eig(Hs*M);
[VHMss,DHMss] = eig(Hs*Ms);
[VMHs,DMHs] = eig(M*Hs);
[VHs,DHs] = eig(Hs);
[VHa,DHa] = eig(Ha);
[VM,DM] = eig(M);
[VMs,DMs] = eig(Ms);
spec_HMs = zeros(Nx,1);
spec_HMss = zeros(Nx,1);
spec_MHs = zeros(Nx,1);
spec_Hs = zeros(Nx,1);
spec_Ha = zeros(Nx,1);
spec_M = zeros(Nx,1);
spec_Ms = zeros(Nx,1);
for n=1:Nx
    spec_HMs(n) = real(DHMs(n,n));
    spec_HMss(n) = real(DHMss(n,n));
    spec_MHs(n) = real(DMHs(n,n));
    spec_Hs(n) = real(DHs(n,n));
    spec_Ha(n) = real(DHa(n,n));
    spec_M(n) = real(DM(n,n));
    spec_Ms(n) = real(DMs(n,n));
end
[spec_HMs,IHMs] = sort(spec_HMs,'descend');
[spec_HMss,IHMss] = sort(spec_HMss,'descend');
[spec_MHs,IMHs] = sort(spec_MHs,'descend');
[spec_Hs,IHs] = sort(spec_Hs,'descend');
[spec_Ha,IHa] = sort(spec_Ha,'descend');
[spec_M,IM] = sort(spec_M);
[spec_Ms,IMs] = sort(spec_Ms);
VHMs(:,:) = VHMs(:,IHMs);
VHMss(:,:) = VHMss(:,IHMss);
VMHs(:,:) = VMHs(:,IMHs);
VHs(:,:) = VHs(:,IHs);
VHa(:,:) = VHa(:,IHa);
VM(:,:) = VM(:,IM);
VMs(:,:) = VMs(:,IMs);


%%% prefactor
mu = abs(spec_HMs(1)+spec_HMs(2));
mu_s = abs(spec_HMss(1)+spec_HMss(2));
lambda = abs(spec_Hs(1)+spec_Hs(2));
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
Epsis = 0.34*10.^linspace(0,log(1e-3)/log(10),10);
epsis = [Epsis(1:9),0.0006,0.0005,0.0004];
% epsis = [0.34,0.2437,0.1747,0.1253,0.0898,0.0644,0.0461,0.0331,0.0237,0.02];


tau_Kramar_mob_var_mass= prefactor_mob_var_mass.*exp(deltaU./epsis);
tau_Kramar_mob_mass= prefactor_mass.*exp(deltaU./epsis);
tau_Kramar_mob = prefactor_mob.*exp(deltaU./epsis);
tau_Kramar_original = prefactor_original.*exp(deltaU./epsis);

epsis_MC = epsis(1:12);
tau_MC = [6.4042e-04,0.0014,0.0027,0.0051,0.0095,0.0161,0.0362,0.0888,0.5623,1.0985,2.4626,9.1712];
% epsis_MC = [epsis(1:14)];
% tau_MC = [7.8672e-04,0.0013,0.0019,0.0028,0.0039,0.0055,0.0079,0.0128,0.0236,0.033,0.0522,0.079];
% tau_MC = [7.9049e-04,0.0013,0.0017,0.0026,0.0037,0.0055,0.0079,0.0129,0.0245,0.036,0.0489,0.0737,0.1832,0.4484];

% dt=1e-7
% tau_MC_1 = [6.9699e-04,0.0011,0.0016,0.0025,0.0035,0.0056,0.0091,0.015,0.0246,0.0381,0.0796,0.1742];

% dt=1e-8
% tau_MC_2 = [7.8672e-04,0.0013,0.0019,0.0028,0.0039,0.0055,0.0079,0.0128,0.0236,0.033,0.0522,0.079,0.1814];

%% plot
semilogy(1./epsis, tau_Kramar_original,"Color", [0.3010 0.7450 0.9330],"LineStyle","--"); hold on
semilogy(1./epsis, tau_Kramar_mob,"k-"); hold on
semilogy(1./epsis, tau_Kramar_mob_mass,'r-'); hold on
semilogy(1./epsis, tau_Kramar_mob_var_mass,'g-'); hold on
semilogy(1./epsis_MC, tau_MC,"k."); hold on
% semilogy(1./[Epsis,0.015,0.012], tau_MC_1,"r*");
% semilogy(1./[Epsis,0.017,0.015,0.012], tau_MC_2,"k*");
title('Eyring-Kramars, h0=1.05', 'interpreter', 'latex','fontsize',16)
legend({'original Eyring-Kramars',...
		'with mobility correction', ...
		'with conserved quantity',...
        'with non-constant mobility',...
        'Monte-Carlo simulation const-mob dt=1e-9'}, 'location',...
	   'best', 'interpreter', 'latex', 'fontsize', 16)
ylabel('mean first passage time $w_B(x\_)$', 'interpreter', 'latex','FontSize',20)
xlabel('noise strength $1/\varepsilon$', 'interpreter', 'latex',"FontSize",20)
grid

% now, save the data for python plot script
save_data_theory = zeros(size(epsis,2),4);
save_data_theory(:,1) = epsis;
save_data_theory(:,2) = tau_Kramar_original;
save_data_theory(:,3) = tau_Kramar_mob;
save_data_theory(:,4) = tau_Kramar_mob_mass;
writematrix(save_data_theory,"STFE_epsis_tau_theory.txt",'Delimiter','space')
save_data_MC = zeros(size(epsis_MC,2),2);
save_data_MC(:,1) = epsis_MC;
save_data_MC(:,2) = tau_MC;
writematrix(save_data_MC,"STFE_epsis_tau_MC.txt",'Delimiter','space')






