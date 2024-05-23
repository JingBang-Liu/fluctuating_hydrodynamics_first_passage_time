clear

load("rp_process_h1.01_full_more.mat")

steps_record_rupture = 5e5;

load("s_temp_1.01.mat")
s = circshift(s,find(s==max(s)),2);

% Spatial and fourier grids
L = 1;
Nx = 128;    % number of nodes
x = linspace(0,L-L/Nx,Nx);
dx = x(2)-x(1);
dk = 2*pi/L;
k = [0:Nx/2,-Nx/2+1:-1]*dk;
epsilon = 0.0005;
U = @(h) sum(1/2.*ifft(1i*k.*fft(h),'symmetric').^2-2/3*pi^2./h.^2)/epsilon*dx;
U_fun = @(h) (1/2.*ifft(1i*k.*fft(h),'symmetric').^2-2/3*pi^2./h.^2)./epsilon;
a = 0*x+mean(s);
dt = 1*1/((8*dk)^4);

U_all = zeros(steps_record_rupture,1);
loss = zeros(steps_record_rupture,1);
for i=1:steps_record_rupture
    U_all(i) = U(full_rp_process(:,i)');
    loss(i) = norm(full_rp_process(:,i)-s');
end
Us = U(s);
ind_loss = find(loss==min(loss));
s_sim = full_rp_process(:,ind_loss)';
disp(ind_loss)

figure(1)
end_cut = 10;
plot(U_all(1:end-end_cut)-(U(s_sim)-Us)-U(a))
hold on
plot(zeros(steps_record_rupture-end_cut,1)+Us-U(a),LineWidth=1.5)
% hold on
% plot(zeros(steps_record_rupture-end_cut,1)+U(s_sim)-U(a))
hold on
plot(zeros(steps_record_rupture-end_cut,1)+U(a)-U(a))

% rp_process_var = squeeze(rp_process_var);
% max_var = max(rp_process_var,[],1);
% hold on
% plot_var = max_var*0.1-6.453;
% plot(plot_var)
% hold on
% plot(plot_var.*0+mean(plot_var(1:9e4)))
ind = find(U_all(1:end-end_cut)>Us,1);
hold on
plot(ones(10,1).*ind_loss,linspace(-1,10,10))
plot(ones(10,1).*3e5,linspace(-1,10,10))
plot(ones(10,1).*(ind_loss+30000),linspace(-1,10,10))
plot(ind_loss+35000,U_all(ind_loss+35000)-U(a),"*")

ylim([U(a)+U(a)/10*epsilon-U(a),Us-Us/10*epsilon-U(a)])


% now, save the data for python plot script
U_all_0005 = zeros(steps_record_rupture+4,1);
U_all_0005(1,1) = Us;
U_all_0005(2,1) = U(s_sim);
U_all_0005(3,1) = U(a);
U_all_0005(4,1) = ind_loss;
U_all_0005(5:end,1) = U_all(:,1);
writematrix(U_all_0005,"STFE_U_all_0.0005.txt",'Delimiter','space')

