ta=500;
tb=500;
pulses=s(k).Pulses;
clust_pos=find(s(k).clusters==cluster); %Finds the cluster position for the clusters we are interested in
clust_time=1000*s(k).times(clust_pos); %Finds the times of the interested cluster: turns into ms
rast=figure;
subplot(3,2,[1:4])
[pointsB,positionB, ~,~]=Raster(pulses,tb,0,clust_time);
set(pointsB(:),'Color',[1 0 0])
hold on
[pointsA,positionA, ~,~]=Raster(pulses,0,ta,clust_time);
set(pointsA(:),'Color',[0 0 1])
title(['Population at Intensity ' num2str(s(k).Intensity{cluster}) '%'])
%xlim([-1*(t_period+500) t_period+500])
xlabel(['Column Number in Base_save is ' num2str(inten_pos(n))])
subplot(3,2,5)
if size(positionB)>0
    plot(1:size(s(k).waveforms(clust_pos(positionB),:),2),...
        s(k).waveforms(clust_pos(positionB),:),'Color',[1 0 0])
end
xlabel(s(k).Intensity{cluster})

subplot(3,2,6)
if size(positionA)>0
        plot(1:size(s(k).waveforms(clust_pos(positionA),:),2),...
            s(k).waveforms(clust_pos(positionA),:),'Color',[0 0 1])
end
xlabel(s(k).Stim{cluster})