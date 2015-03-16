%Find the mean of the 

tbef=500;
ta=500
sigDiff=zeros(size(normptsh(:,3+gauss_size:end-gauss_size)));

for k=1:size(normptsh,1)
    meanBase=mean(normptsh(k,1+tbase-tbef:tbase-1))
    stdBase=std(normptsh(k,1+tbase-tbef:tbase-1))./sqrt(length(normptsh(k,1+tbase-tbef:tbase-1)))
    greater=find(normptsh(k,3+gauss_size:end-gauss_size)>meanBase+2*stdBase);
    lesser=find(normptsh(k,3+gauss_size:end-gauss_size)<meanBase-2*stdBase);%...
    %zspot=find(normptsh(k,3+gauss_size:end-gauss_size)==0);
    if length(greater)>0
        sigDiff(k,greater)=1;
    end
    if length(lesser)>0
        sigDiff(k,lesser)=-1;
    end
%     if length(zspot)>0
%         sigDiff(k,zspot)=-.5;
%     end
end

figure
imagesc(sigDiff)
colormap('hot')

%Plot spots where there is significant difference
figure;plot(mean(abs(sigDiff)))
%Shorten matrix to 110ms after the TMS pulse
sigDiffafter=sigDiff(:,500:500+1+110);
%plot shortened matrix
figure;imagesc(abs(sigDiffafter))

%Plot the sorted mean for each row, this will help remove cells without
%effect
figure;plot(sort(mean(abs(sigDiffafter),2)))
%Find all cells where the mean is zero for the time frame
pos=find(mean(abs(sigDiffafter),2)>0);

sigCells=sigDiffafter(pos,:);
blockdata=AreaDate(pos,:);
figure;imagesc(sigCells)
colorbar
figure;plot(mean(sigCells<0))
hold on;plot(mean(sigCells>0),'g')

figure
for n=1:9
    posinten=find(cell2mat(blockdata(:,4))<=n*10 & ...
        cell2mat(blockdata(:,4))>10*(n-1));
    subplot(3,3,n)
    plot(mean(sigCells(posinten,:)<0))
    hold on;plot(mean(sigCells(posinten,:)>0),'g')
    title([num2str(n*10) '% with Cell Count= ' num2str(length(posinten))])
end
    