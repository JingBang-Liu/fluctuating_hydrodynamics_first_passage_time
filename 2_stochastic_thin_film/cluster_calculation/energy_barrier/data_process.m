% combine the files
clear
steps_record_rupture = 5e5;
Nseeds = 40;
Nevents = 10;
Nx = 128;
full_rp_process = zeros(Nx,steps_record_rupture);
rp_process_var = zeros(Nx,steps_record_rupture);
count = zeros(1,steps_record_rupture);
count_disp = 0;
tol = 0.01;

for ll=1:Nseeds
    for kk = 1:Nevents
        fn = "rp_process_h1.01_seed_"+num2str(ll)+"_rpcount_"+num2str(kk)+".mat";
        load(fn);
        mass = mean(rupture_process,1);
        ind = find(mass==0,1,"last");
        if isempty(ind)
            ind = 0;
        end
        count(1,(ind+1):end) = count(1,(ind+1):end) + 1;
        ind2 = find(count==0,1,"last");
        ind3 = max([ind,ind2]);
        count_disp = count_disp + 1;
        disp("count "+num2str(count_disp)+" ind "+num2str(ind)+" ind3 "+num2str(ind3))
        for ii = 1:Nx
            full_rp_process(ii,(ind3+1):end) = (full_rp_process(ii,(ind3+1):end)...
            .*(count(1,(ind3+1):end)-1)+rupture_process(ii,(ind3+1):end))./count(1,(ind3+1):end);
        end
    end
end
count_disp = 0;
ind2 = find(count==0,1,"last");
if isempty(ind2)
    ind2 = 0
end
for ll=1:Nseeds
    for kk=1:Nevents
        fn = "rp_process_h1.01_seed_"+num2str(ll)+"_rpcount_"+num2str(kk)+".mat";
        load(fn);
        count_disp = count_disp + 1;
        disp(count_disp)
        for ii=1:Nx
            rp_process_var(ii,ind2+1:end) = rp_process_var(ii,ind2+1:end)...
            +(rupture_process(ii,ind2+1:end)-full_rp_process(ii,ind2+1:end)).^2./ ...
            count(1,ind2+1:end);
        end
    end
end

save("rp_process_h1.01_full_more.mat","full_rp_process","rp_process_var","count")

exit
