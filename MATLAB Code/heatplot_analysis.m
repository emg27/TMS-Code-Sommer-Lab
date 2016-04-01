clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
close all
clf

brainloc = zeros(size(normptsh,1),1);
for n = 1:size(normptsh,1)
brainloc(n,1) = AreaDate{n,3};
end

reds = linspace(0,1,9)';
% greens = [0:0.25:1 0.75:-0.25:0]';
greens = zeros(9,1);
blues = linspace(1,0,9)';
 cmap = [reds greens blues];

figure
pSh=find(normptsh(:,1)==0 & brainloc(:)==1);
pSt=find(normptsh(:,1)==1& brainloc(:)==1);
shamps=normptsh(pSh,:);
stimps=normptsh(pSt,:);
%shamps=normptsh(normptsh(:,1)==0,:);
%stimps=normptsh(normptsh(:,1)==1,:);
peakpoints=zeros(9,2);
peaktimes=peakpoints;
for n=1:9
pos=find(shamps(:,2)<=n*10 & shamps(:,2)>10*(n-1));
subplot(6,3,n);imagesc(shamps(pos,3+gauss_size:end-gauss_size),[-1 1])
title(['Sham ' num2str(n) '0%'])
line([tbase+1 tbase+1],[0 length(pos)+1],'Color','k')
subplot(6,3,n+9)
avgSham(n,:)=nanmean(shamps(pos,3+gauss_size:end-gauss_size));
avgShamSM(n,:)=smooth(avgSham(n,:),25);
allsham=shamps(pos,3+gauss_size:end-gauss_size);
tempSh=nanmax(allsham(:,tbase:tbase+101)');
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

plot(avgSham(n,:),'r');
hold on
%plot(nanmean(allptsh(pSh(pos),3:end)),'b')
plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
line([tbase+1 tbase+1],[-1 1],'Color','k')
ylim([-.25 .6])
%line([500 500],[-.01 0.03],'Color','k')
title(['Sham ' num2str(n) '0%'])
xlim([1 size(avgSham(n,:),2)])
peakpoints(n,1)=max(avgSham(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
peaktimes(n,1)=find(avgSham(n,:)==peakpoints(n,1))-tbase;
peakpointsSM(n,1)=max(avgShamSM(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
peaktimesSM(n,1)=find(smooth(avgSham(n,:),25)==peakpointsSM(n,1))-tbase;
end
figure
for n=1:9
pos=find(stimps(:,2)<=n*10 & stimps(:,2)>10*(n-1));
subplot(6,3,n);imagesc(stimps(pos,3+gauss_size:end-gauss_size),[-1 1])
title(['Stim ' num2str(n) '0%'])
line([tbase+1 tbase+1],[0 length(pos)+1],'Color','k')
subplot(6,3,n+9)
avgStim(n,:)=nanmean(stimps(pos,3+gauss_size:end-gauss_size));
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

plot(avgStim(n,:),'r');
hold on
%plot(nanmean(allptsh(pSt(pos),3:end)),'b')
plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
line([tbase+1 tbase+1],[-1 1],'Color','k')
ylim([-.25 .6])
%line([500 500],[-.01 0.03],'Color','k')
title(['Stim ' num2str(n) '0%'])
xlim([1 size(avgSham(n,:),2)])
%xlim([15 1008])
peakpoints(n,2)=max(avgStim(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
%temp=findpeaks(avgStim(n,3+tbase+1:tbase+101));
%peakpoints(n,2)=max(temp);
peakpointsSM(n,2)=max(avgStimSM(n,tbase+1:tbase+101)); %Finds and stores the max point after the TMS pulse
peaktimes(n,2)=find(avgStim(n,:)==peakpoints(n,2))-tbase;
peaktimesSM(n,2)=find(smooth(avgStim(n,:),25)==peakpointsSM(n,2))-tbase;
end

figure
subplot(2,1,1)
plot(peakpoints(1:9,:),'o'); legend('Sham','Stim')
title(['Raw Peak points, Gauss=' num2str(gauss_size) 'ms'])
subplot(2,1,2)
plot(peakpointsSM(1:9,:),'o'); legend('Sham','Stim')
title('Smoothed by 25 points')

figure
subplot(2,1,1)
plot(peaktimes(1:9,:),'o'); legend('Sham','Stim') %Raw Plot times
title(['Raw Times, Gauss=' num2str(gauss_size) 'ms'])
subplot(2,1,2)
plot(peaktimesSM(1:9,:),'o'); legend('Sham','Stim') %Raw Plot times
title('Smooth Peak Points by 25 25 points')

figure
subplot(2,1,1)
for p = 1:9
plot(-50:ta,avgSham(p,(tbase+1)-50:(tbase+1)+ta)','Color',cmap(p,:));hold on
end
line([0 0],[-.2 .5],'Color','k')
%legend('10','20','30','40','50','60','70','80','90')
title(['Sham, Gauss=' num2str(gauss_size) 'ms'])
xlim([-50 200])
%axis([tbase-43 tbase+ta+3 -.15 .4])
subplot(2,1,2)
for q = 1:9
plot(-50:ta,avgStim(q,(tbase+1)-50:(tbase+1)+ta)','Color',cmap(q,:)); hold on
end
line([0 0],[-.2 .5],'Color','k')
legend('10','20','30','40','50','60','70','80','90')
title(['Stim, Gauss=' num2str(gauss_size) 'ms'])
%axis([tbase-43 tbase+ta+3 -.15 .4])
xlim([-50 200])

figure
subplot(1,5,[1 2])
imagesc(allptsh(pSt,3+gauss_size:end-gauss_size))
line([tbase+1 tbase+1],[0 length(pSt)+1],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
title('Raw Stim')
xlim([0 tbase+ta])
colorbar
subplot(1,5,[3 4])
imagesc(normptsh(pSt,3+gauss_size:end-gauss_size))
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
imagesc(allptsh(pSh,3+gauss_size:end-gauss_size))
line([tbase+1 tbase+1],[0 length(pSh)+1],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
title('Raw Sham')
xlim([0 tbase+ta])
colorbar
subplot(1,5,[3 4])
imagesc(normptsh(pSh,3+gauss_size:end-gauss_size))
line([tbase+1 tbase+1],[0 length(pSh)+1],'Color','k')
ylabel('Cells-All Intensities')
title('Normalized Sham')
xlabel('Time (ms)')
%title('Norm Sham')
colorbar
subplot(1,5,5)
imagesc(cell2mat(AreaDate(pSh,3)))
