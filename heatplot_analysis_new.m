clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
[colormapfile, pathnamecolor]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathnamecolor colormapfile])
close all

figure
pSh=find(normptsh(:,1)==0);
pSt=find(normptsh(:,1)==1);
shamps=normptsh(pSh,:);
stimps=normptsh(pSt,:);
%shamps=normptsh(normptsh(:,1)==0,:);
%stimps=normptsh(normptsh(:,1)==1,:);
peakpoints=zeros(9,2);
peaktimes=peakpoints;

CmapHeatPlotcheck=CmapHeatPlot6; %Use either 4 or 6
Tbefpulse=200;
time=-(3+gauss_size+tbase):(tbase+gauss_size-1);
plotrange=[3+gauss_size+tbase+1-Tbefpulse:size(shamps,2)+1-gauss_size];
for n=1:9
pos=find(shamps(:,2)<=n*10 & shamps(:,2)>10*(n-1));

%subplot(6,3,n);imagesc(shamps(pos,3+gauss_size:end-gauss_size),[-1 1])
subplot(6,3,n);imagesc(shamps(pos,plotrange),[-1 1])
colormap(CmapHeatPlotcheck) %jet

title(['Sham ' num2str(n) '0%'])
line([Tbefpulse+1 Tbefpulse+1],[0 length(pos)+1],'Color','k')
subplot(6,3,n+9)
avgSham(n,:)=nanmean(shamps(pos,:));
avgShamSM(n,:)=smooth(avgSham(n,:),25);
allsham=shamps(pos,3+gauss_size:end-gauss_size);
tempSh=nanmax(allsham(:,tbase:tbase+101)');
%colormap jet
for j=1:length(pos)
    if mean(allsham(j,tbase:tbase+101)-tempSh(j))~=0
        posMaxSh=find(allsham(j,tbase:tbase+101)==tempSh(j));
        maxAllSh(j)=posMaxSh(1);
    else
        %figure;plot(allsham(j,tbase:tbase+101)')
        maxAllSh(j)=nan;
    end
end
temp2Sh(n)=nanmean(maxAllSh);

plot(time(plotrange),avgSham(n,plotrange),'r');
hold on
%plot(nanmean(allptsh(pSh(pos),3:end)),'b')
plot(time(plotrange),0*(time(plotrange)),'k--')
%line([Tbefpulse+1 Tbefpulse+1],[-1 1],'Color','k')
line([0 0],[-1 1],'Color','k')
ylim([-.25 .6])
title(['Sham ' num2str(n) '0%'])
%xlim([1 size(avgSham(n,plotrange),2)])
xlim([time(plotrange(1)) time(plotrange(end))])
peakpoints(n,1)=max(avgSham(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
peaktimes(n,1)=find(avgSham(n,:)==peakpoints(n,1))-tbase;
peakpointsSM(n,1)=max(avgShamSM(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
peaktimesSM(n,1)=find(smooth(avgSham(n,:),25)==peakpointsSM(n,1))-tbase;
end
figure
for n=1:9
pos=find(stimps(:,2)<=n*10 & stimps(:,2)>10*(n-1));
subplot(6,3,n);imagesc(stimps(pos,plotrange),[-1 1])
colormap(CmapHeatPlotcheck)
title(['Stim ' num2str(n) '0%'])
line([Tbefpulse+1 Tbefpulse+1],[0 length(pos)+1],'Color','k')
subplot(6,3,n+9)
avgStim(n,:)=nanmean(stimps(pos,:));
avgStimSM(n,:)=smooth(avgStim(n,:),25);
allstim=stimps(pos,3+gauss_size:end-gauss_size);
tempSt=nanmax(allstim(:,tbase:tbase+101)');
for j=1:length(pos)
    if mean(allstim(j,tbase:tbase+101)-tempSt(j))~=0
        posMaxSt=find(allstim(j,tbase:tbase+101)==tempSt(j));
        maxAllSt(j)=posMaxSt(1);
    else
        %figure;plot(allstim(j,tbase:tbase+101)')
        maxAllSt(j)=nan;
    end
end
temp2St(n)=nanmean(maxAllSt);

plot(time(plotrange),avgStim(n,plotrange),'r');
hold on
%plot(nanmean(allptsh(pSt(pos),3:end)),'b')
plot(time(plotrange),0*(time(plotrange)),'k--')
line([0 0],[-1 1],'Color','k')
ylim([-.25 .6])
%line([500 500],[-.01 0.03],'Color','k')
title(['Stim ' num2str(n) '0%'])
xlim([time(plotrange(1)) time(plotrange(end))])
%xlim([15 1008])
peakpoints(n,2)=max(avgStim(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
%temp=findpeaks(avgStim(n,3+tbase+1:tbase+101));
%peakpoints(n,2)=max(temp);
peakpointsSM(n,2)=max(avgStimSM(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
peaktimes(n,2)=find(avgStim(n,:)==peakpoints(n,2))-tbase;
peaktimesSM(n,2)=find(smooth(avgStim(n,:),25)==peakpointsSM(n,2))-tbase;
end

% figure
% subplot(2,1,1)
% plot(peakpoints(1:9,:),'o'); legend('Sham','Stim')
% title(['Raw Peak points, Gauss=' num2str(gauss_size) 'ms'])
% subplot(2,1,2)
% plot(peakpointsSM(1:9,:),'o'); legend('Sham','Stim')
% title('Smoothed by 25 points')
% 
% figure
% subplot(2,1,1)
% plot(peaktimes(1:9,:),'o'); legend('Sham','Stim') %Raw Plot times
% title(['Raw Times, Gauss=' num2str(gauss_size) 'ms'])
% subplot(2,1,2)
% plot(peaktimesSM(1:9,:),'o'); legend('Sham','Stim') %Raw Plot times
% title('Smooth Peak Points by 25 25 points')

figure
subplot(2,1,1)
plot(-50:ta,avgSham(1:9,(tbase+1)-50:(tbase+1)+ta)','LineWidth',2);
hold on
line([0 0],[-.2 .5],'Color','k')
%legend('10','20','30','40','50','60','70','80','90')
title(['Sham, Gauss=' num2str(gauss_size) 'ms'])
xlim([-50 200])
%axis([tbase-43 tbase+ta+3 -.15 .4])
subplot(2,1,2)
plot(-50:ta,avgStim(1:9,(tbase+1)-50:(tbase+1)+ta)','LineWidth',2);
line([0 0],[-.2 .5],'Color','k')
legend('10','20','30','40','50','60','70','80','90')
title(['Stim, Gauss=' num2str(gauss_size) 'ms'])
%axis([tbase-43 tbase+ta+3 -.15 .4])
xlim([-50 200])

avgStim_sub = mean(avgStim(1:5,:));
std_Stim_sub = std(avgStim(1:5,:));
avgStim_supra = mean(avgStim(6:9,:));
std_Stim_supra = std(avgStim(6:9,:));
avgSham_sub = mean(avgSham(1:5,:));
std_Sham_sub = std(avgSham(1:5,:));
avgSham_supra = mean(avgSham(6:9,:));
std_Sham_supra = std(avgSham(6:9,:));
figure
subplot(2,1,1)
crop = (tbase+1)-50:(tbase+1)+ta;
plot(-50:ta,avgSham_sub(crop),'b')
%plot_variance(-50:ta,avgSham_sub(crop)-std_Sham_sub(crop),avgSham_sub(crop)+std_Sham_sub(crop),'b');
hold on
plot(-50:ta,avgSham_supra(crop),'r');
%plot_variance(-50:ta,avgSham_supra(crop)-std_Sham_supra(crop),avgSham_supra(crop)+std_Sham_supra(crop),'r');
line([0 0],[-.2 .5],'Color','k')
title(['Sham, Gauss=' num2str(gauss_size) 'ms'])
xlim([-50 200])
subplot(2,1,2)
plot(-50:ta,avgStim_sub(crop),'b')
%plot_variance(-50:ta,avgStim_sub(crop)-std_Stim_sub(crop),avgStim_sub(crop)+std_Stim_sub(crop),'b');
hold on
plot(-50:ta,avgStim_supra(crop),'r');
%plot_variance(-50:ta,avgStim_supra(crop)-std_Stim_supra(crop),avgStim_supra(crop)+std_Stim_supra(crop),'r');
line([0 0],[-.2 .5],'Color','k')
legend('Subthreshold','Suprathreshold')
title(['Stim, Gauss=' num2str(gauss_size) 'ms'])
xlim([-50 200])

figure
subplot(1,5,1:2)
imagesc(allptsh(pSt,3+gauss_size:end-gauss_size)); colormap jet
line([tbase+1 tbase+1],[0 length(pSt)+1],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
title('Raw Stim')
xlim([0 tbase+ta])
colorbar
subplot(1,5,[3 4])
imagesc(normptsh(pSt,3+gauss_size:end-gauss_size)); colormap jet
line([tbase+1 tbase+1],[0 length(pSt)+1],'Color','k')
ylabel('Cells-All Intensities')
title('Normalized Stim')
xlabel('Time (ms)')
%title('Norm Sham')
colorbar
subplot(1,5,5)
imagesc(cell2mat(AreaDate(pSt,3)))

figure
subplot(1,5,1:2)
imagesc(allptsh(pSh,3+gauss_size:end-gauss_size)); colormap jet
line([tbase+1 tbase+1],[0 length(pSh)+1],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
title('Raw Sham')
xlim([0 tbase+ta])
colorbar
subplot(1,5,[3 4])
imagesc(normptsh(pSh,3+gauss_size:end-gauss_size)); colormap jet
line([tbase+1 tbase+1],[0 length(pSh)+1],'Color','k')
ylabel('Cells-All Intensities')
title('Normalized Sham')
xlabel('Time (ms)')
%title('Norm Sham')
colorbar
subplot(1,5,5)
imagesc(cell2mat(AreaDate(pSh,3)))
