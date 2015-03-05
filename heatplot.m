clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
close all

tbase=500; %amount of baseline collected
ta=500; %amount of time after the TMS pulse
gauss_size=2.5;

counter=0; %Will count the number of cells used 
allptsh=[]; 
normptsh=[];
%2012 Data: 1 to 113
%2013 Data: 114 to 287
%2014 Data: 288 to end
for k=1:size(s,2)%1:113%size(s,2)
    if(length(s(k).Pulses)>0) & median(diff(s(k).Pulses))>4 %& ...
            %size(s(k).Stim,1)>0 & strcmp(s(k).Stim(1),'Stim')==1 &...
            %size(s(k).Intensity,1)>0 & strcmp(s(k).Intensity(1),'90')==1
        pulses=s(k).Pulses;
        firerate=s(k).FireRate;
        for g=1:max(s(k).clusters)
            cluster=find(s(k).clusters==g);
            if length(cluster)>0 & length(cluster)>length(pulses)
                figure(10)
                [spk_d,trl_fr,bin_start_times,baseline,mean_trl_fr]=...
                    psth1block(pulses,tbase+gauss_size,ta+gauss_size, 1000*s(k).times(cluster), gauss_size,0);
                close(10);
                if size(s(k).Stim,1)>0 & strcmp(s(k).Stim(1),'Stim')==1
                    stim=1;
                elseif size(s(k).Stim,1)>0 & strcmp(s(k).Stim(1),'Sham')==1
                    stim=0;
                else
                    stim=3;
                end
                if size(s(k).Intensity,1)>0
                    inten=s(k).Intensity(1);
                else
                    inten=-10;   
                end
                if size(inten)==0
                    continue
                end
                if size(nonzeros(mean_trl_fr),1)<1 | strcmp(s(k).BrainArea,'V1') | strcmp(s(k).BrainArea,'VI')
                    continue
                end
                counter=counter+1;
                AreaDate{counter,1}=s(k).BrainArea;
                AreaDate{counter,2}=s(k).Name;
                if strcmp(s(k).BrainArea,'FEF')
                    AreaDate{counter,3}=-1;
                elseif strcmp(s(k).BrainArea,'V1') | strcmp(s(k).BrainArea,'VI')
                    AreaDate{counter,3}=-.5;
                elseif strcmp(s(k).BrainArea,'M1') | strcmp(s(k).BrainArea,'MI')
                    AreaDate{counter,3}=0;
                else
                    AreaDate{counter,3}=1;
                end
                AreaDate{counter,4}=inten;
%                 stim=1;
%                 inten=nan;
                allptsh=[allptsh; stim inten mean_trl_fr];
                normptsh=[normptsh; stim inten mean_trl_fr/max(abs(mean_trl_fr))];
            end
        end
    end
end
figure
subplot(1,7,[1 2])
imagesc(allptsh(:,3+gauss_size:end-gauss_size))
line([tbase tbase],[0 counter+1],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
xlim([0 tbase+ta])
%title('Stim')
% figure
subplot(1,7,[3 4])
imagesc(normptsh(:,3+gauss_size:end-gauss_size))
line([tbase tbase],[0 counter+1],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
%title('Norm Sham')
subplot(1,7,5)
imagesc(cell2mat(AreaDate(:,3)))
subplot(1,7,6)
imagesc(normptsh(:,1))
subplot(1,7,7)
imagesc(normptsh(:,2))
[newfilename,saveloc]=uiputfile('*.mat'); %Change to the location where you want to save to
save([saveloc newfilename],'normptsh','allptsh','AreaDate',...
    'tbase','ta','gauss_size')


% figure
% pSh=find(normptsh(:,1)==0);
% pSt=find(normptsh(:,1)==1);
% shamps=normptsh(pSh,:);
% stimps=normptsh(pSt,:);
% %shamps=normptsh(normptsh(:,1)==0,:);
% %stimps=normptsh(normptsh(:,1)==1,:);
% peakpoints=zeros(10,2);
% peaktimes=peakpoints;
% for n=1:10
% pos=find(shamps(:,2)<=n*10 & shamps(:,2)>10*(n-1));
% subplot(4,5,n);imagesc(shamps(pos,3:end),[-1 1])
% title(['Sham ' num2str(n) '0%'])
% line([500 500],[0 length(pos)+1],'Color','k')
% subplot(4,5,n+10)
% avgSham(n,:)=nanmean(shamps(pos,3:end));
% avgShamSM(n,:)=smooth(avgSham(n,:),25);
% plot(avgSham(n,:),'r');
% hold on
% %plot(nanmean(allptsh(pSh(pos),3:end)),'b')
% plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
% line([500 500],[-1 1],'Color','k')
% ylim([-.25 .6])
% %line([500 500],[-.01 0.03],'Color','k')
% title(['Sham ' num2str(n) '0%'])
% xlim([15 1008])
% peakpoints(n,1)=max(avgSham(n,tbase+3:tbase+203)); %Finds and stores the max point after the TMS pulse
% peaktimes(n,1)=find(avgSham(n,:)==peakpoints(n,1))-tbase;
% peakpointsSM(n,1)=max(smooth(avgSham(n,tbase+3:tbase+203),25)); %Finds and stores the max point after the TMS pulse
% peaktimesSM(n,1)=find(smooth(avgSham(n,:),25)==peakpointsSM(n,1))-tbase;
% end
% figure
% for n=1:10
% pos=find(stimps(:,2)<=n*10 & stimps(:,2)>10*(n-1));
% subplot(4,5,n);imagesc(stimps(pos,3:end),[-1 1])
% title(['Stim ' num2str(n) '0%'])
% line([500 500],[0 length(pos)+1],'Color','k')
% subplot(4,5,n+10)
% avgStim(n,:)=nanmean(stimps(pos,3:end));
% avgStimSM(n,:)=smooth(avgStim(n,:),25);
% plot(avgStim(n,:),'r');
% hold on
% %plot(nanmean(allptsh(pSt(pos),3:end)),'b')
% plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
% line([500 500],[-1 1],'Color','k')
% ylim([-.25 .6])
% %line([500 500],[-.01 0.03],'Color','k')
% title(['Stim ' num2str(n) '0%'])
% xlim([15 1008])
% peakpoints(n,2)=max(avgStim(n,tbase+3:tbase+203)); %Finds and stores the max point after the TMS pulse
% peakpointsSM(n,2)=max(smooth(avgStim(n,tbase+3:tbase+203),25)); %Finds and stores the max point after the TMS pulse
% peaktimes(n,2)=find(avgStim(n,:)==peakpoints(n,2))-tbase;
% peaktimesSM(n,2)=find(smooth(avgStim(n,:),25)==peakpointsSM(n,2))-tbase;
% end
% 
% figure
% subplot(2,1,1)
% plot(peakpoints(1:9,:),'o'); legend('Sham','Stim')
% title(['Raw Peak points, Gauss=' num2str(gauss_size) 'ms'])
% subplot(2,1,2)
% plot(peakpointsSM(1:9,:),'o'); legend('Sham','Stim')
% title('Smoothed by 25 points')
% figure
% subplot(2,1,1)
% plot(peaktimes(1:9,:),'o'); legend('Sham','Stim') %Raw Plot times
% title(['Raw Times, Gauss=' num2str(gauss_size) 'ms'])
% subplot(2,1,2)
% plot(peaktimesSM(1:9,:),'o'); legend('Sham','Stim') %Raw Plot times
% title('Smooth Peak Points by 25 25 points')
% figure
% subplot(2,1,1)
% plot(avgSham(1:9,:)');hold on
% line([tbase+3 tbase+3],[-.15 .4],'Color','k')
% %legend('10','20','30','40','50','60','70','80','90')
% title(['Sham, Gauss=' num2str(gauss_size) 'ms'])
% xlim([tbase-43 tbase+ta+3])
% %axis([tbase-43 tbase+ta+3 -.15 .4])
% subplot(2,1,2)
% plot(avgStim(1:9,:)')
% line([tbase+3 tbase+3],[-.15 .4],'Color','k')
% legend('10','20','30','40','50','60','70','80','90')
% title(['Stim, Gauss=' num2str(gauss_size) 'ms'])
% %axis([tbase-43 tbase+ta+3 -.15 .4])
% xlim([tbase-43 tbase+ta+3])