%%Calculate the baselines for cells before and after the tms pulse and then
%%plot the neuron response, after vs before.

%load data
file='spike_base_stim_10.mat';
load(file);

%set the parameter
t_period=200; %Time period for average firing rate before and after
base_bef=[];
base_aft=[];
figure
%Run through the structure
for n=1:size(spike_base,2)
    blockdata=s(spike_base(1,n)); %save the particular block we are interested in
    pulses=blockdata.Pulses;
    clust_pos=find(blockdata.clusters==spike_base(2,n)); %Finds the cluster position for the clusters we are interested in
    clust_time=1000*blockdata.times(clust_pos); %Finds the times of the interested cluster: turns into ms
    [~,~, num_spike_bef,~]=Raster(pulses,t_period,0,clust_time);
    [~,~, num_spike_aft,~]=Raster(pulses,0,t_period,clust_time);
    close
    fire_bef=num_spike_bef./(size(pulses,1)*t_period/1000);
    fire_aft=num_spike_aft./(size(pulses,1)*t_period/1000);
    base_bef=[base_bef sum(fire_bef)];
    base_aft=[base_aft sum(fire_aft)];
end
plot(base_bef,base_aft,'o')
hold on
unitx=linspace(0,ceil(1.1*max(base_bef)),1000);
plot(unitx,unitx,'k-')
xlim([0 ceil(1.1*max(base_bef))])