% combine the files
clear
full_events = zeros(96,128);
full_waiting_time = zeros(96,1);

for ll=1:8
    fn = "h01.01_0.0004_"+num2str(ll)+".mat";
    load(fn);
    full_waiting_time(((ll-1)*12+1):ll*12,1) = waiting_time(:,1);
    full_events(((ll-1)*12+1):ll*12,:) = events(:,:);
end

save("h01.01_0.0004_full.mat","full_events","full_waiting_time")