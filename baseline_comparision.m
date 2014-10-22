%%Calculate the baselines for cells before and after the tms pulse and then
%%plot the neuron response, after vs before.
close all
%load data
file='base_save.mat';
load(file);

%set the parameter
t_period=200; %Time period for average firing rate before and after
base_bef=[];
base_aft=[];

%Create a file for each intensity
intensity=unique(base_save(3,:));
for h=1:size(intensity,2)
    inten_pos=find(base_save(3,:)==intensity(h));
    pop(h,:)=size(inten_pos);
    base_inten=base_save(:,inten_pos);
    base_bef=nan(25,pop(h,2));
    base_aft=nan(25,pop(h,2));
    %Run through the structure
    for n=1:size(base_inten,2)
        blockdata=s(base_inten(1,n)); %save the particular block we are interested in
        pulses=blockdata.Pulses;
        clust_pos=find(blockdata.clusters==base_inten(2,n)); %Finds the cluster position for the clusters we are interested in
        clust_time=1000*blockdata.times(clust_pos); %Finds the times of the interested cluster: turns into ms
        figure(99)
        [~,~, num_spike_bef,~]=Raster(pulses,t_period,0,clust_time);
        [~,~, num_spike_aft,~]=Raster(pulses,0,t_period,clust_time);
        close(99)
        fire_bef=num_spike_bef./(t_period/1000);
        fire_aft=num_spike_aft./(t_period/1000);
%         base_bef=[base_bef mean(fire_bef(1))];
%         base_aft=[base_aft mean(fire_aft(1))];
        base_bef(1:length(fire_bef),n)=fire_bef;
        base_aft(1:length(fire_bef),n)=fire_aft;
    end
    figure
    plot(base_bef(1,:),base_aft(1,:),'o')
    hold on
    unitx=linspace(min(min(base_bef)),ceil(1.1*max(max(base_bef))),1000);
    plot(unitx,unitx,'k-')
    %xlim([-ceil(1.1*max(log10(base_bef)))  ceil(1.1*max(log10(base_bef)))])
    xlabel(['Scale Baseline ' num2str(t_period) ' Before TMS Pulse']);
    ylabel(['Baseline ' num2str(t_period) ' After TMS Pulse']);
    title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(pop(h,2))]);
    figure
    hist(base_aft(1,:)-base_bef(1,:))
    title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(pop(h,2))]);
end