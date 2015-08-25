function[spikes,time,spikesM,spikesMV]=centerspks(spks,firerate,full,check)
%%center all the waveforms in the system

baseline=mean(mean(spks(:,1:15))); %Finds DC offset
spks=spks-baseline; %Remove DC offset
minspks=min(spks(:,ceil(size(spks,2)./3):end-5),[],2);

time=nan(1,full);
posmin=find(mean(spks)==min(mean(spks)));
time(full/2-posmin:full/2-posmin-1+size(spks,2))=1000*([1:size(spks,2)]-posmin)./firerate;
spikes=[];

if check==1

for n=1:size(spks,1)
    spkstemp=nan(1,full);
   % spikesV=nan(1,300);

    pos=find(spks(n,:)==minspks(n));
    spkstemp(full/2-pos(1):full/2-pos(1)-1+size(spks,2))=spks(n,:);
    
    spikes=[spikes;spkstemp];
end

spikesM=nanmean(spikes);
spikesM(isnan(spikesM))=0;
spikesMV=nanstd(spikes);
spikesMV=spikesMV;
spikesMV(isnan(spikesMV))=0;
else
    spikesM=nan(1,full);
    spikesMV=nan(1,full);
    spikesM(full/2-posmin:full/2-posmin-1+size(spks,2))=mean(spks);
    spikesMV(full/2-posmin:full/2-posmin-1+size(spks,2))=std(spks);
end


