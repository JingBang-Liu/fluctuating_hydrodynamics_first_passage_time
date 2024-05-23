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
deltaU = (Us - Ua)*0.83;
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
mu = abs(spec_HMs(1));
mu_s = abs(spec_HMss(1));
lambda = abs(spec_Hs(1));
nhat = VM(:,1);
nhat_s = VMs(:,1);
nHsn = abs(nhat'*inv(Hs)*nhat);
nHan = abs(nhat'*inv(Ha)*nhat);
nHsn_s = abs(nhat_s'*inv(Hs)*nhat_s);
nHan_s = abs(nhat_s'*inv(Ha)*nhat_s);
spec_Hs_new = spec_Hs;
spec_Hs_new(2) = 1;
prefactor_mob_var_mass = pi/mu_s*sqrt(abs(prod(spec_Hs_new./spec_Ha)))*sqrt(nHsn_s/nHan_s); % conserved quantity and variable mobility
prefactor_mass = pi/mu*sqrt(abs(prod(spec_Hs_new./spec_Ha)))*sqrt(nHsn/nHan); % with conserved quantity
prefactor_mob = pi/mu*sqrt(abs(prod(spec_Hs_new./spec_Ha))); % with mobility
prefactor_original = pi/lambda*sqrt(abs(prod(spec_Hs_new./spec_Ha))); % original
disp("prefactor original "+num2str(prefactor_original))
disp("prefactor mobility "+num2str(prefactor_mob))
disp("prefactor mass     "+num2str(prefactor_mass))
disp("prefactor mob var  "+num2str(prefactor_mob_var_mass))

 %% calculate time
Epsis = 0.34*10.^linspace(0,log(1e-3)/log(10),10);
epsis = [Epsis(1:9),0.0006,0.0005,0.0004];


tau_Kramar_mob_var_mass= prefactor_mob_var_mass.*exp(deltaU./epsis);
tau_Kramar_mob_mass= prefactor_mass.*exp(deltaU./epsis);
tau_Kramar_mob = prefactor_mob.*exp(deltaU./epsis);
tau_Kramar_original = prefactor_original.*exp(deltaU./epsis);

epsis_MC = epsis(1:12);
tau_MC = [6.3967e-04,0.0014,0.0027,0.005,0.0093,0.0152,0.0322,0.0689,0.3257,0.5640,1.1065,3.3739];

%% plot
semilogy(1./epsis, tau_Kramar_original,"Color", [0.3010 0.7450 0.9330],"LineStyle","--"); hold on
semilogy(1./epsis, tau_Kramar_mob,"Color",[0 0.4470 0.7410],"LineStyle","--"); hold on
semilogy(1./epsis, tau_Kramar_mob_mass,'r-'); hold on
semilogy(1./epsis, tau_Kramar_mob_var_mass,'g-'); hold on
semilogy(1./epsis_MC, tau_MC,"k."); hold on
title('Eyring-Kramars, h0=1.01', 'interpreter', 'latex','fontsize',16)
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


