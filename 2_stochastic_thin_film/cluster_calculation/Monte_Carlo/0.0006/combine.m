% combine the files
clear
full_events = zeros(100,128);
full_waiting_time = zeros(100,1);

for ll=1:4
    fn = "h01.01_0.0006_"+num2str(ll)+".mat";
    load(fn);
    full_waiting_time(((ll-1)*25+1):ll*25,1) = waiting_time(:,1);
    full_events(((ll-1)*25+1):ll*25,:) = events(:,:);
end

save("h01.01_0.0006_full.mat","full_events","full_waiting_time")