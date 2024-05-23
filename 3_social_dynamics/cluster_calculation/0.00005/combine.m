% combine the files
clear
full_events = zeros(100,64);
full_waiting_time = zeros(100,1);

for ll=1:4
    fn = "rho0_0.79076_eps_5e-05_"+num2str(ll)+".mat";
    load(fn);
    full_waiting_time(((ll-1)*25+1):ll*25,1) = waiting_time(:);
    full_events(((ll-1)*25+1):ll*25,:) = events(:,:);
end

save("eps_0.00005.mat","full_events","full_waiting_time")
