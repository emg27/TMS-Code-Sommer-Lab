%%Calculate the baselines for cells before and after the tms pulse and then
%%plot the neuron response, after vs before.
close all
%load data
file='base_save_SFN_stim_rtms.mat';
load(file);

%set the parameter
t_period=200; %Time period for average firing rate before and after
base_bef=[];
base_aft=[];

dose=figure;
%Create a file for each intensity
intensity=unique(base_save(3,:));
for h=size(intensity,2)-1
    inten_pos=find(base_save(3,:)==intensity(h));
    pop(h,:)=size(inten_pos);
    base_inten=base_save(:,inten_pos);
    base_bef=nan(70,pop(h,2));
    base_aft=nan(70,pop(h,2));
    lastpulse=[];
    freq=[];
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
        [~,fr_bef,bin_bef,~,~]=psth1block(pulses,t_period+60,60,clust_time,30,1);
        pos_bef=find(bin_bef>=-t_period & bin_bef<0);
        fire_bef=mean(fr_bef(:,pos_bef),2);
        [~,fr_aft,bin_aft,~,~]=psth1block(pulses,60,t_period+60,clust_time,30,1);
        pos_aft=find(bin_aft>=0 & bin_aft<t_period);
        fire_aft=mean(fr_aft(:,pos_aft),2);
%         [~,~, num_spike_bef,~]=Raster(pulses,t_period,0,clust_time);
%         hold on
%         [~,~, num_spike_aft,~]=Raster(pulses,0,t_period,clust_time);
         close(rast)
         
%          figure;
%         subplot(3,2,[1:4])
%         [pointsB,positionB, ~,~]=Raster(pulses(1),t_period+500,0,clust_time);
%         set(pointsB(:),'Color',[1 0 0])
%         hold on
%         [pointsA,positionA, ~,~]=Raster(pulses(end),0,t_period+500,clust_time);
%         set(pointsA(:),'Color',[0 0 1])
%         title(['Population at Intensity ' num2str(intensity(h)) '%, Cell=' ...
%             num2str(n) ' of ' num2str(size(base_bef,2))])
%         xlabel(['Column Number in Base_save is ' num2str(inten_pos(n))])  
%         xlim([-1*(t_period+500) t_period+500])
%         %ylim([0 2])
%         subplot(3,2,5)
%         if size(positionB)>0
%         plot(1:size(blockdata.waveforms(clust_pos(positionB),:),2),...
%             blockdata.waveforms(clust_pos(positionB),:),'Color',[1 0 0])
%         end
%         subplot(3,2,6)
%         if size(positionA)>0
%         plot(1:size(blockdata.waveforms(clust_pos(positionA),:),2),...
%             blockdata.waveforms(clust_pos(positionA),:),'Color',[0 0 1])
%         end
%         fire_bef=num_spike_bef./(t_period/1000);
%         fire_aft=num_spike_aft./(t_period/1000);
%         base_bef(1:length(fire_bef),n)=fire_bef;
%         base_aft(1:length(fire_bef),n)=fire_aft;
        base_bef(1:length(fire_bef),n)=fire_bef*1000;
        base_aft(1:length(fire_bef),n)=fire_aft*1000;
        lastpos=find(diff(pulses)>2000);
        pulse_diff=diff(pulses);
        if size(lastpos,1)<=0
            lastpulse=[lastpulse fire_aft(end)*1000];
            freq=[freq 1000./mean(pulse_diff)];
        else
            lastpulse=[lastpulse fire_aft(lastpos(1))*1000];
            freq=[freq 1000./mean(pulse_diff(1:lastpos-1))];
        end
        
        
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
    base_all=[lastpulse base_bef(1,:)];
    unitx=linspace(min(base_all),ceil(1.1*max(base_all)),1000);
    %unitx=linspace(0,40,1000);
    tmsmark=[];
    for k=1:length(lastpulse)
        if freq(k)>=9 & freq(k)<=11
            plot(base_bef(1,k),lastpulse(k),'o','MarkerFaceColor',[1 0 0])
            tmsmark(k)=1;
        elseif freq(k)>=4 & freq(k)<=6
            plot(base_bef(1,k),lastpulse(k),'o','MarkerFaceColor',[0 1 0])
            tmsmark(k)=2;
        elseif freq(k)>=.4 & freq(k)<=1.6
            plot(base_bef(1,k),lastpulse(k),'o','MarkerFaceColor',[0 0 1])
            tmsmark(k)=3;
        else
            tmsmark(k)=0;
        end
    end
    plot(unitx,unitx,'k-')
    %xlim([min(unitx) max(unitx)])
    %ylim([min(unitx) max(unitx)])
    %xlim([-ceil(1.1*max(log10(base_bef)))  ceil(1.1*max(log10(base_bef)))])
    xlabel(['Scale Baseline ' num2str(t_period) 'ms Before TMS Pulse']);
    ylabel(['Baseline ' num2str(t_period) 'ms After TMS Pulse']);
    title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(size(base_bef(1,:),2))]);
    base_change=lastpulse-base_bef(1,:);

%     figure
%     plot(base_change,'go','MarkerFaceColor',[0 1 0])
%     hold on
%     plot(base_aft(1,:),'bo','MarkerFaceColor',[0 0 1])
%     plot(base_bef(1,:),'ro','MarkerFaceColor',[1 0 0])
%     legend('Difference','After','Before',0)
%     title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(base_bef(1,base_bef(1,:)>0)))]);


    pulse=1;
    base_change=base_aft(pulse,:)-base_bef(pulse,:);
    base_last=zeros(1,size(base_aft,2));
    %find the baseline after the last pulse
    for cells=1:size(base_aft,2)
        posLP=isfinite(base_aft(:,cells));
        base_last(cells)=base_aft(find(posLP,1,'last'),cells);
    end
    
    figure
    subplot(1,2,1)
    edge=linspace(0,15,21);
    hist(base_change,15);%,edge);
    [Hchange, pchange]=ttest(base_aft(pulse,:));
    %[pchange, Hchange]=ranksum(base_change,normpdf(length(base_change(pulse,:))));
    %axis([min(edge) max(edge) 0 20])
%     ylim([0 45])
    title(['Pulse: ' num2str(pulse) ' Population at Intensity ' num2str(intensity(h)) '% Stim, Cell Count=' num2str(size(base_bef,2))]);
    xlabel(['Effect Comparision=' num2str(Hchange) ' p-value is ' num2str(pchange)])
%     xlabel(['% Change After/Before Cell Count=' num2str(size(base_bef(pulse,base_bef(pulse,:)<=0),2))])
    subplot(1,2,2)
%     hist(base_last-base_bef(pulse,:),20);
    hist(base_bef(pulse,:),15);
    [Hbef, pvbef]=ttest(base_bef(pulse,:));
%     [Hbef, pvbef]=ttest(base_last-base_bef(pulse,:));%,normpdf(length(base_bef(pulse,:))));
%[pvbef, Hbef]=ranksum(base_bef(pulse,:),normpdf(length(base_bef(pulse,:))));
%     ylim([0 60])
    %set(f,'FaceColor','g');
%     axis([min(edge2) max(edge2) 0 30])
     xlabel(['Baseline Last-First Comparision=' num2str(Hbef) ' p-value is ' num2str(pvbef)])
    
    pulse=1;
    figure
%     subplot(1,2,1)
    %hist(base_change,linspace(min(base_change),max(base_change),10))
    hist(lastpulse-base_bef(pulse,:),20)
%     edge=linspace(0,15,21);
%     %[N,Bin]=histc(base_change./(base_bef(1,:)+.0001),edge);
%     [N,Bin]=histc(lastpulse(base_bef(pulse,:)>0)./(base_bef(pulse,base_bef(pulse,:)>0)),edge);
%     bar(edge,N,'histc')
%     axis([min(edge) max(edge) 0 15])
    title(['Pulse: ' num2str(pulse) ' Population at Intensity ' num2str(intensity(h)) '% Stim, Cell Count=' num2str(size(base_bef,2))]);
    xlabel(['% Change After/Before Cell Count=' num2str(size(base_bef(1,base_bef(1,:)>0),2))])
%     subplot(1,2,2)
%     %edge2=linspace(0,15,21);
%     [N2,Bin2]=histc(lastpulse(base_bef(pulse,:)<=0),edge2);
%     if size(N2,2)>0
%         f=bar(edge2,N2,'histc')
%         set(f,'FaceColor','g');
%         axis([min(edge2) max(edge2) 0 10])
%         xlabel(['Spontaneous Firing Cell Count=' num2str(size(base_bef(pulse,base_bef(pulse,:)<=0),2))])
%     end
    
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
% figure
% plot(base_bef(:,1),'o')
% title(['Population at Intensity ' num2str(intensity(h)) '%, Cell Count=' num2str(length(base_bef(1,base_bef(1,:)>0)))]);

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
    %all_cell_regress_bef=[all_cell_regress_bef [baseregressbef; length(num_pulse_bef)]];
    [rho_bef,p_bef]=corr(num_pulse_bef',base_pulse_bef,'type','Spearman')
    
    X_diff=[];
    X_diff=[num_pulse_diff' ones(size(base_pulse_bef,1),1)];
    baseregressdiff=regress(base_pulse_diff,X_diff);
    
    [rho_diff,p_diff]=corr(num_pulse_diff',base_pulse_diff,'type','Spearman')
    
   % if neuron==26%h==5 & (neuron==9 | neuron==44 | neuron==40) % p_diff<=0.05
%         figure
%         subplot(2,1,1)
%         plot(base_pulse_diff,'o')
%         title(['Population at Intensity ' num2str(intensity(h)) '%, Cell=' ...
%             num2str(neuron) ' of ' num2str(size(base_bef,2))]);
%         hold on
%         plot(num_pulse_diff,X_diff*baseregressdiff,'r-')
%         subplot(2,1,2)
%         plot(cell2mat(wave(neuron,2)),cell2mat(wave(neuron,1)))
%         title(['P-value= ' num2str(p_diff)])
   % end
    if p_diff<=0.05
        all_cell_regress_diff=[all_cell_regress_diff [baseregressdiff; length(num_pulse_diff); 1]];
    else
        all_cell_regress_diff=[all_cell_regress_diff [baseregressdiff; length(num_pulse_diff); 0]];
    end
    
    if p_bef<=0.05
        all_cell_regress_bef=[all_cell_regress_bef [baseregressbef; length(num_pulse_bef); 1]];
    else
        all_cell_regress_bef=[all_cell_regress_bef [baseregressbef; length(num_pulse_bef); 0]];
    end
end

baser=figure
baseB=figure
hist_cell_slopes=nan(size(all_cell_regress_bef,2),3);
hist_cell_base_slopes=nan(size(all_cell_regress_bef,2),3);
for k=1:size(all_cell_regress_bef,2)
% subplot(2,1,1)
% plot(k,all_cell_regress_bef(1,k),'o',...
%     'Color',[0 0 all_cell_regress_bef(4,k)],...
%     'MarkerFaceColor',[all_cell_regress_bef(4,k) 0 0])
%  hold on
% title(sprintf('Population Regress Slope for the Baseline %d ms before a TMS pulse at Intensity %d', t_period,intensity(h)))
% xlabel(['Cell Count= ' num2str(size(base_bef,2))])
% %ylim([-2.5 4])
if tmsmark(k)>0
figure(baser)
subplot(3,1,tmsmark(k))
% plot(k,all_cell_regress_bef(1,k).*(all_cell_regress_bef(3,k)-1),'o',...
%     'Color',[0 0 all_cell_regress_bef(4,k)],...
%     'MarkerFaceColor',[all_cell_regress_bef(4,k) 0 0])
 hold on
% subplot(2,2,3)
plot(k,all_cell_regress_diff(1,k),'o',...
    'Color',[0 0 all_cell_regress_diff(4,k)],...
    'MarkerEdgeColor',(tmsmark(k)==1).*[1 0 0]+(tmsmark(k)==2).*[0 1 0]+...
    (tmsmark(k)==3).*[0 0 1],...
    'MarkerFaceColor',[0 0 all_cell_regress_diff(4,k)])
title('Difference Between before and after')

figure(baseB)
subplot(3,1,tmsmark(k))
% plot(k,all_cell_regress_bef(1,k).*(all_cell_regress_bef(3,k)-1),'o',...
%     'Color',[0 0 all_cell_regress_bef(4,k)],...
%     'MarkerFaceColor',[all_cell_regress_bef(4,k) 0 0])
 hold on
% subplot(2,2,3)
plot(k,all_cell_regress_bef(1,k),'o',...
    'Color',[0 0 all_cell_regress_bef(4,k)],...
    'MarkerEdgeColor',(tmsmark(k)==1).*[1 0 0]+(tmsmark(k)==2).*[0 1 0]+...
    (tmsmark(k)==3).*[0 0 1],...
    'MarkerFaceColor',[0 0 all_cell_regress_bef(4,k)])
title('Baseline Changes')
%  hold on
% 
%     
% subplot(2,2,4)
% plot(k,all_cell_regress_diff(1,k).*(all_cell_regress_diff(3,k)-1),'o',...
%     'Color',[0 0 all_cell_regress_diff(4,k)],...
%     'MarkerFaceColor',[0 0 all_cell_regress_diff(4,k)])
%  hold on
hist_cell_slopes(k,tmsmark(k))=all_cell_regress_diff(1,k);
hist_cell_base_slopes(k,tmsmark(k))=all_cell_regress_bef(1,k);
end
end
figure(baseB)
subplot(3,1,1)
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
ylim([-0.6 0.6])
subplot(3,1,2)
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
ylim([-0.6 0.6])
subplot(3,1,3)
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
ylim([-0.6 0.6])

figure(baser)
subplot(3,1,1)
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
ylim([-0.2 0.2])
%ylim([-0.6 0.65])
subplot(3,1,2)
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
ylim([-0.2 0.2])
%ylim([-0.6 0.65])
subplot(3,1,3)
plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
ylim([-0.2 0.2])
%ylim([-0.6 0.65])
% subplot(2,2,3)
% rgslope=all_cell_regress_diff(1,:);
% plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
% xlabel(sprintf('Regress <0: %d and Regress >0: %d Regress=0: %d',...
%     length(rgslope(rgslope<0)), length(rgslope(rgslope>0)),...
%     length(rgslope(rgslope==0))))
% subplot(2,2,4)
% plot(0:size(base_bef,2)+1,zeros(size(base_bef,2)+2,1),'g-')
%hist(rgslope,20)

figure
% subplot(2,1,1)
% hist(all_cell_regress_bef(1,:),15)
% title('Slope of Baselines before pulse')
% xlabel(['Intensity= ' num2str(intensity(h))])
% subplot(2,1,2)
subplot(1,3,1)
hist(hist_cell_slopes(:,1),linspace(-.2,.2,20))
title('Difference Between before and after')
%xlim([-0.2 0.2])
 axis([-.2 0.2 0 18])
xlabel(['Intensity= ' num2str(intensity(h))])
subplot(1,3,2)
hist(hist_cell_slopes(:,2),linspace(-.2,.2,20))
xlim([-0.2 0.2])
 axis([-.2 0.2 0 18])
%axis([-.2 0.2 0 12])
%axis([-.1 0.6 0 20])
title('Difference Between before and after')
subplot(1,3,3)
hist(hist_cell_slopes(:,3),linspace(-.2,.2,20))
xlim([-0.2 0.2])
 axis([-.2 0.2 0 18])
%axis([-.2 0.2 0 12])
%axis([-.1 0.6 0 20])
title('Difference Between before and after')


figure
% subplot(2,1,1)
% hist(all_cell_regress_bef(1,:),15)
% title('Slope of Baselines before pulse')
% xlabel(['Intensity= ' num2str(intensity(h))])
% subplot(2,1,2)
subplot(1,3,1)
hist(hist_cell_base_slopes(:,1),linspace(-.6,.6,20))
xlim([-0.6 0.6])
ylim([0 25])
%axis([-.2 0.2 0 18])
title('Baseline Frequency: 10' )
xlabel(['Intensity= ' num2str(intensity(h))])
subplot(1,3,2)
hist(hist_cell_base_slopes(:,2),linspace(-.6,.6,20))
xlim([-0.6 0.6])
ylim([0 25])
%axis([-.2 0.2 0 18])
title('Baseline Frequency: 5' )
subplot(1,3,3)
hist(hist_cell_base_slopes(:,3),linspace(-.6,.6,20))
xlim([-0.6 0.6])
ylim([0 25])
%axis([-.2 0.2 0 18])
title('Baseline Frequency: 1' )

figure(dose)
subplot(1,2,1)
zscore=(all_cell_regress_bef(1,:)-mean(all_cell_regress_bef(1,:)))./std(all_cell_regress_bef(1,:));
plot(intensity(h),zscore,'o')
%plot(intensity(h),mean(all_cell_regress_bef(1,:)),'o')
%errorbar(intensity(h),mean(all_cell_regress_bef(1,:)),var(all_cell_regress_bef(1,:)),'o')
hold on
xlim([0 100])
subplot(1,2,2)
plot(intensity(h),mean(all_cell_regress_diff(1,:)),'o')
errorbar(intensity(h),mean(all_cell_regress_diff(1,:)),var(all_cell_regress_diff(1,:)),'o')
hold on
xlim([0 100])
figure
hist(zscore)
title(['Z-score for Before BAseline regression slope, Intensity=' num2str(intensity(h))])
end
figure(dose)
subplot(1,2,1)
title('Baseline')
xlabel('Intensity')
ylabel('Slope value')
subplot(1,2,2)
title('Effect')
xlabel('Intensity')
ylabel('Slope value')