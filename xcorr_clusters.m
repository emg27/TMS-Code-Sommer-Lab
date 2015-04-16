%close all
intensity=90; 
Tstart=500;
%gauss_size=5;
%pos=find(normptsh(:,1)<2);% & normptsh(:,2)>=60);% & normptsh(:,2)>intensity-10 & normptsh(:,2)<=intensity);
pos=find(normptsh(:,1)==1 & normptsh(:,2)>=70);
words='Stim greater than or equal to 60';
N=length(pos);
cutoff=1.6;

corrMat=zeros(N);

waves=normptsh(pos,3+gauss_size+Tstart+1:end-gauss_size-400);
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
xlabel(words)

joy=figure
joy2=figure
joy3=figure
k1=0;
for k=1:max(groups)
posB=find(groups==k);
if length(posB)<2
    continue
end
k1=k1+1;
figure(joy2)
subplot(2,3,k1)
StShInten=normptsh(pos(posB),1:2);
Stimpts=StShInten(StShInten(:,1)==1);
Shampts=StShInten(StShInten(:,1)==0);
bar([1 0],[length(Stimpts)/sum(normptsh(pos,1)==1),...
    length(Shampts)/sum(normptsh(pos,1)==0)])
xlabel(words)
figure(joy3)
subplot(2,3,k1)
hist(normptsh(pos(posB),2)+100*normptsh(pos(posB),1),6)
xlabel(words)
%Plot the waveforms
files=cell2mat(AreaDate(pos(posB),6));
neuron=cell2mat(AreaDate(pos(posB),5));
%figure
%spikepl=[];
% for n=1:length(posB);
%     spikes=find(s(files(n)).clusters==neuron(n));
%     waveforms=s(files(n)).waveforms(spikes,:);
%     meanWave=mean(waveforms);
%     check(n)=length(meanWave);
%     %time=1000*(0:length(meanWave)-1)/s(files(n)).FireRate;
%     [spikes,time,spikesM,spikesMV]=centerspks(waveforms,s(files(n)).FireRate,300,0);
%     if n<=25
%         n1=n;
%     elseif n==26
%         figure
%         n1=n-25;
%     else
%         n1=n-25;
%     end
%     subplot(5,5,n1)
%     %plot(spikesM/max(abs(spikesM)))
%     plot(time,spikesM/max(abs(spikesM)))
%     hold on
%     spikepl=[spikepl; spikesM];
%     
%     title(['Group ' num2str(k) 'Neuron ' num2str(n)])
%     xlabel([num2str(files(n)) ' ' AreaDate{pos(posB(n)),2}])
%     ylim([-1.1 1.1])
%     xlim([-1.5 3])
% end
% % figure(joy2)
% % subplot(2,5,k)
% % plot(nanmean(spikepl))
% % title(['Group ' num2str(k) 'Neuron ' num2str(n)])
% %     xlabel([num2str(files(n)) ' ' AreaDate{pos(posB(n)),2}])
% %     ylim([-1.1 1.1])
% %     xlim([-1.5 3])
figure(joy)
subplot(2,3,k1)
% figure
hold on
imagesc(normptsh(pos(posB),3+gauss_size:end-gauss_size),[-1 1]);
line([501 501],[0 length(posB)+1],'Color','k')
title(['Group ' num2str(k)])
axis([0 1000 0 length(posB)+1])
xlabel(words)
end

