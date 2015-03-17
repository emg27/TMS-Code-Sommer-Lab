close all
intensity=80;
Tstart=500;
%gauss_size=5;
pos=find(normptsh(:,2)>intensity-10 & normptsh(:,2)<=intensity & normptsh(:,1)==1);
N=length(pos);
cutoff=1.32;

corrMat=zeros(N);

waves=normptsh(pos,3+gauss_size+Tstart+1:end-gauss_size);
%waves=normptsh(pos,3+gauss_size:end-gauss_size);
%waves=sigDiff(pos,Tstart+1:end);
temp=[];
for p=1:N
    for q=1:N
        temp=xcorr(waves(p,:),waves(q,:),'coeff');
        %covar2(p,q)=max(temp);
        %maxval=temp(find(abs(temp)==max(abs(temp))));
        maxval=temp(ceil(size(temp,2)/2))+0.1E-8;
        corrMat(p,q)=maxval(1);
    end
end

corrMat = tril(corrMat,-1);
dissimilarity = 1 - corrMat(find(corrMat))';
Z = linkage(dissimilarity,'complete');
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
%groups = cluster(Z,'cutoff',cutoff);
figure
dendrogram(Z,0,'colorthreshold',cutoff)
joy=figure
%joy2=figure
for k=1:max(groups)
posB=find(groups==k);

%Plot the waveforms
files=cell2mat(AreaDate(pos(posB),6));
neuron=cell2mat(AreaDate(pos(posB),5));
figure
spikepl=[];
for n=1:length(posB);
    spikes=find(s(files(n)).clusters==neuron(n));
    waveforms=s(files(n)).waveforms(spikes,:);
    meanWave=mean(waveforms);
    check(n)=length(meanWave);
    %time=1000*(0:length(meanWave)-1)/s(files(n)).FireRate;
    [spikes,time,spikesM,spikesMV]=centerspks(waveforms,s(files(n)).FireRate,300,0);
    subplot(6,4,n)
    %plot(spikesM/max(abs(spikesM)))
    plot(time,spikesM/max(abs(spikesM)))
    hold on
    spikepl=[spikepl; spikesM];
    
    title(['Group ' num2str(k) 'Neuron ' num2str(n)])
    xlabel([num2str(files(n)) ' ' AreaDate{pos(posB(n)),2}])
    ylim([-1.1 1.1])
    xlim([-1.5 3])
end
% figure(joy2)
% subplot(2,5,k)
% plot(nanmean(spikepl))
% title(['Group ' num2str(k) 'Neuron ' num2str(n)])
%     xlabel([num2str(files(n)) ' ' AreaDate{pos(posB(n)),2}])
%     ylim([-1.1 1.1])
%     xlim([-1.5 3])
figure(joy)
subplot(3,3,k)
% figure
hold on
imagesc(normptsh(pos(posB),3+gauss_size:end-gauss_size),[-1 1]);
line([501 501],[0 length(posB)+1],'Color','k')
title(['Group ' num2str(k)])
axis([0 1000 0 length(posB)+1])
end
