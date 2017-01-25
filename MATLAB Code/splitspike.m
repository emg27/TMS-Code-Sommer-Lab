%% Prepare Space


%% Variables
s_proxy = s;
totInd = 1;


%% Analysis
for i = 1:length(s)
    clusterProxy = s(i).clusters;
    totClusters = unique(clusterProxy(find(clusterProxy ~= 0 & clusterProxy ~= 9)));
    
    for j = 1:length(totClusters)
        clusterInd = find(clusterProxy == totClusters(j));
        
        s_proxy(totInd).Name = s(i).Name;
        s_proxy(totInd).Blackout = s(i).Blackout;
        s_proxy(totInd).Animal = s(i).Animal;
        s_proxy(totInd).FireRate = s(i).FireRate;
        s_proxy(totInd).Pulses = s(i).Pulses;
        s_proxy(totInd).times = s(i).times(clusterInd);
        s_proxy(totInd).clusters = s(i).clusters(clusterInd);
        s_proxy(totInd).waveforms = s(i).waveforms(clusterInd,:);
        s_proxy(totInd).Stim = s(i).Stim;
        s_proxy(totInd).Intensity = s(i).Intensity;
        s_proxy(totInd).Area = s(i).Area;
        s_proxy(totInd).Depth = s(i).Depth;
        s_proxy(totInd).BrainArea = s(i).BrainArea;
        s_proxy(totInd).Locate = s(i).Locate;
        s_proxy(totInd).ID = totClusters(j);
        totInd = totInd + 1;
    end
end

s = s_proxy;