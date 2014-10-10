function[spikes,spikesV,spikesM,spikesMV]=centerspks(spks,full)
%%center all the waveforms in the system

minspks=min(spks(:,30:90),[],2);

spikes=[];

for n=1:size(spks,1)
    spkstemp=nan(1,full);
   % spikesV=nan(1,300);

    pos=find(spks(n,:)==minspks(n));
    spkstemp(full/2-pos:full/2-pos-1+size(spks,2))=spks(n,:);
    
    spikes=[spikes;spkstemp];
end

spikesM=nanmean(spikes);
spikesM(isnan(spikesM))=0;
spikesV=nanstd(spikes);
spikesMV=spikesV;
spikesMV(isnan(spikesV))=0;


