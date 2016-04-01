% K means clustering
% Number of clusters is chosen before generating/viewing clusters
% 
% clear
% [filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
% load([pathname filename])
% close all
clf

Tstart = 500;
Tend = Tstart+150; 

brainloc = zeros(size(normptsh,1),1);
blackout = zeros(size(normptsh,1),1);
animal = zeros(size(normptsh,1),1);
for n = 1:size(normptsh,1)
    brainloc(n,1) = AreaDate{n,3};
    blackout(n,1) = str2double(AreaDate{n,8});
    if strcmp(AreaDate{n,7},'Magneto')==1
        animal(n,1) = 1;
    elseif strcmp(AreaDate{n,7},'Tacoma')==1
        animal(n,1) = -1;
    end
end

% Identify neurons of interest
% responsive = find_responsive_blocks(normptsh,gauss_size); % could be allptsh or normptsh
neurons = find(normptsh(:,1)==1 & brainloc(:)==1 & blackout(:)<=3); % & responsive(:)==0); % & animal(:)==-1);
words = 'M1 Sham, Blackout<=3ms';

% Isolate activity profiles of interest, in this case, first 150 ms after TMS pulse
waves=normptsh(neurons,3+gauss_size+1+Tstart:3+gauss_size+1+Tend);

% K means clustering
% Evaluate range of K values to find optimal K
% Optimal is defined by the chosen criterion (CalinskiHarabasz, DaviesBouldin, silhouette, gap)
kmax = 10;
reps = 100;
clust = zeros(size(waves,1),kmax);
for k=1:kmax
    rng('default')
    clust(:,k) = kmeans(waves,k,'Start','sample','MaxIter',100,'Replicates',reps);
    % Start: chooses initial cluster centroid positions
        % sample: chooses initial centroid positions at random
        % plus: chooses initial centroid positions based on kmeans++ algorithm
    % MaxIter: iterates clustering using those centroid positions
        % changes the clustering MaxIter times until it reaches a global minimum
        % (min sum of distances between points and their respective centroids)
    % Replicates: re-chooses new set of initial cluster centroid positions
        % refers to completely separate runs of the clustering algorithm
end
eva = evalclusters(waves,clust,'DaviesBouldin');
numgroups = eva.OptimalK;
idx = clust(:,numgroups);

counter=0;
pop_psth = [];
percent_inten_cluster = zeros(numgroups,9);

% Plot results
for n=1:numgroups
pos = find(idx==n);
if length(pos)<2
    continue
end
counter=counter+1;

figure(1)
subplot(3,7,counter)
hold on
imagesc(normptsh(neurons(pos),3+gauss_size:end-gauss_size),[-1 1]);
line([501 501],[0 length(pos)+1],'Color','k')
title(['Group ' num2str(n)])
axis([0 1000 0 length(pos)+1])
xlabel(words)

pop_psth(n,:) = nanmean(normptsh(neurons(pos),3+gauss_size:end-gauss_size));
subplot(3,7,counter+7)
hold on
plot(pop_psth(n,:),'r')
plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
line([tbase+1 tbase+1],[-1 1],'Color','k')
axis([0 1000 -0.5 1])
title(['Group ' num2str(n)]) 
xlabel(words)

for j = 1:9
    inten = length(find(normptsh(neurons,2)==j*10));
        % calculates number of blocks at each intensity
    inten_cluster = length(find(normptsh(neurons(pos),2)==j*10));
        % calculates number of blocks in this cluster at each intensity
    percent_inten_cluster(counter,j) = 100*inten_cluster./inten;
end

subplot(3,7,counter+2*7)
bar(10:10:90,percent_inten_cluster(counter,:))
title(['Group ' num2str(n)])
axis([0 100 0 100])
ylabel('Percentage') % percent of all blocks at this intensity that assign to this cluster
xlabel(words)

end

