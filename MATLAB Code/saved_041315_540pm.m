maxPeaks=max(normptsh(:,3+gauss_size+tbase+1:3+gauss_size+tbase+102)');

for n=1:size(maxPeaks,2)
    if mean(normptsh(n,3+gauss_size+tbase+1:3+gauss_size+tbase+102)-maxPeaks(n))~=0
        peaktimes=find(normptsh(n,3+gauss_size+tbase+1:3+gauss_size+tbase+102)==maxPeaks(n));
        timePeaks(n)=peaktimes(1);
    else
        timePeaks(n)=nan;
    end
end