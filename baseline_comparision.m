%%Calculate the baselines for cells before and after the tms pulse and then
%%plot the neuron response, after vs before.
close all
%load data
file='base_save_SFN.mat';
load(file);

%set the parameter
t_period=300; %Time period for average firing rate before and after
base_bef=[];
base_aft=[];

dose=figure;
%Create a file for each intensity
intensity=unique(base_save(3,:));
for h=1:9;%size(intensity,2)
    inten_pos=find(base_save(3,:)==intensity(h));
    pop(h,:)=size(inten_pos);
    base_inten=base_save(:,inten_pos);
    base_bef=nan(25,pop(h,2));
    base_aft=nan(25,pop(h,2));
    wave=wave_save(inten_pos,:);
    %Run through the structure
    for n=1:size(base_inten,2)
        blockdata=s(base_inten(1,n)); %save the particular block we are interested in
        pulses=blockdata.Pulses;
        clust_pos=find(blockdata.clusters==base_inten(2,n)); %Finds the cluster position for the clusters we are interested in
        clust_time=1000*blockdata.times(clust_pos); %Finds the times of the interested cluster: turns into ms
        rast=figure;
        subplot(2,1,1)
        title([blockdata.Name ' Intensity: ' num2str(h*10)])
        [~,~, num_spike_bef,~]=Raster(pulses,t_period,0,clust_time);
        hold on
        [~,~, num_spike_aft,~]=Raster(pulses,0,t_period,clust_time);
        close(rast)
        fire_bef=num_spike_bef./(t_period/1000);
        fire_aft=num_spike_aft./(t_period/1000);
        base_bef(1:length(fire_bef),n)=fire_bef;
        base_aft(1:length(fire_bef),n)=fire_aft;
%         if base_aft(1,n)-base_bef(1,n)>50
%             subplot(2,1,2)
%             plot(cell2mat(wave(n,2)),cell2mat(wave(n,1)))
%             xlabel(sprintf('Fire Before: %d Fire After: %d Diff: %d',...
%                 base_bef(1,n), base_aft(1,n), base_aft(1,n)-base_bef(1,n)))
%         else
%             close(rast)
%         end
    end
    figure
    
    hold on
    base_all=[base_bef base_aft];
    %unitx=linspace(min(min(base_all(1,:))),ceil(1.1*max(max(base_all(1,:)))),1000);
    unitx=linspace(0,125,1000);
    plot(base_bef(1,:),base_aft(1,:),'o','MarkerFaceColor',[0 0 1])
    plot(unitx,unitx,'k-')
    %xlim([min(unitx) max(unitx)])
    %ylim([min(unitx) max(unitx)])
    %xlim([-ceil(1.1*max(log10(base_bef)))  ceil(1.1*max(log10(base_bef)))])
    xlabel(['Scale Baseline ' num2str(t_period) 'ms Before TMS Pulse']);
    ylabel(['Baseline ' num2str(t_period) 'ms After TMS Pulse']);
    title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(size(base_bef,2))]);
    base_change=base_aft(1,:)-base_bef(1,:);

%     figure
%     plot(base_change,'go','MarkerFaceColor',[0 1 0])
%     hold on
%     plot(base_aft(1,:),'bo','MarkerFaceColor',[0 0 1])
%     plot(base_bef(1,:),'ro','MarkerFaceColor',[1 0 0])
%     legend('Difference','After','Before',0)
%     title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(base_bef(1,base_bef(1,:)>0)))]);
    pulse=1;
    figure
    subplot(1,2,1)
    %hist(base_change,linspace(min(base_change),max(base_change),10))
    edge=linspace(0,15,21);
    %[N,Bin]=histc(base_change./(base_bef(1,:)+.0001),edge);
    [N,Bin]=histc(base_aft(pulse,base_bef(pulse,:)>0)./(base_bef(pulse,base_bef(pulse,:)>0)),edge);
    bar(edge,N,'histc')
    axis([min(edge) max(edge) 0 20])
    title(['Pulse: ' num2str(pulse) ' Population at Intensity ' num2str(intensity(h)) '% Stim, Cell Count=' num2str(size(base_bef,2))]);
    xlabel(['% Change After/Before Cell Count=' num2str(size(base_bef(pulse,base_bef(pulse,:)<=0),2))])
    subplot(1,2,2)
    edge2=linspace(0,40,21);
    [N2,Bin2]=histc(base_aft(pulse,base_bef(pulse,:)<=0),edge2);
    f=bar(edge2,N2,'histc');
    set(f,'FaceColor','g');
    axis([min(edge2) max(edge2) 0 30])
    xlabel(['Spontaneous Firing Cell Count=' num2str(size(base_bef(pulse,base_bef(pulse,:)<=0),2))])
    
%     figure
%     subplot(3,1,1)
%     plot(1:25,nanmean(base_aft-base_bef,2),'o-')
%     hold on
%     plot(1:25,nanmedian(base_aft-base_bef,2),'ro-')
%     axis([0 20 -8 8])
%     title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(pop(h,2))]);
%     subplot(3,1,2)
%     plot(1:25,base_aft-base_bef,'o')
%     hold on
%     plot(1:25,zeros(size(1:25)),'k-')
%     xlim([0 20])
%     subplot(3,1,3)
%     plot(1:24,diff(base_aft-base_bef),'o')
%     xlim([0 20])
all_cell_regress_bef=[];
all_cell_regress_diff=[];
for neuron=1:size(base_bef,2)
    base_pulse_bef=base_bef(:,neuron);
    base_pulse_bef=base_pulse_bef(isfinite(base_pulse_bef));
    num_pulse_bef=1:size(base_pulse_bef,1);
    
    base_pulse_diff=base_aft(:,neuron)-base_bef(:,neuron);
    base_pulse_diff=base_pulse_diff(isfinite(base_pulse_diff));
    num_pulse_diff=1:size(base_pulse_diff,1);
    
    X_bef=[];
    X_bef=[num_pulse_bef' ones(size(base_pulse_bef,1),1)];
    baseregressbef=regress(base_pulse_bef, X_bef);
    all_cell_regress_bef=[all_cell_regress_bef baseregressbef];
    
    X_diff=[];
    X_diff=[num_pulse_diff' ones(size(base_pulse_bef,1),1)];
    baseregressdiff=regress(base_pulse_diff,X_diff);
    all_cell_regress_diff=[all_cell_regress_diff baseregressdiff];
    
%     figure
%     plot(base_pulse_diff,'o')
%     title(['Population at Intensity ' num2str(intensity(h)) '%, Cell=' ...
%         num2str(neuron) ' of ' num2str(size(base_bef,2))]);
%     hold on
%     plot(num_pulse_diff,X_diff*baseregressdiff,'r-')
end

figure
subplot(3,1,1)
plot(all_cell_regress_bef(1,:),'o')
hold on
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'k-')
title(sprintf('Population Regress Slope for the Baseline %d ms before a TMS pulse at Intensity %d', t_period,intensity(h)))
xlabel(['Cell Count= ' num2str(size(base_bef,2))])
%ylim([-2.5 4])

subplot(3,1,2)
plot(all_cell_regress_diff(1,:),'o')
title('Difference Between before and after')
rgslope=all_cell_regress_diff(1,:);
xlabel(sprintf('Regress <0: %d and Regress >0: %d Regress=0: %d',...
    length(rgslope(rgslope<0)), length(rgslope(rgslope>0)),...
    length(rgslope(rgslope==0))))
subplot(3,1,3)
hist(rgslope,20)

figure(dose)
subplot(1,2,1)
plot(intensity(h),mean(all_cell_regress_bef(1,:)),'o')
errorbar(intensity(h),mean(all_cell_regress_bef(1,:)),var(all_cell_regress_bef(1,:)),'o')
hold on
xlim([0 100])
subplot(1,2,2)
plot(intensity(h),mean(all_cell_regress_diff(1,:)),'o')
errorbar(intensity(h),mean(all_cell_regress_diff(1,:)),var(all_cell_regress_diff(1,:)),'o')
hold on
xlim([0 100])
end