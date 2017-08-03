% K means clustering
% Number of clusters is chosen before generating/viewing clusters

% First run heatplot code
    % heatplot code receives struct 's' of data, generates 'normptsh'
% Now run this k means code
    % k means code receives 'normptsh', generates k means figures

clf

Tstart = 500;
Tend = Tstart+200; % 200 ms after TMS pulse

brainloc = zeros(size(normptsh,1),1);
blackout = zeros(size(normptsh,1),1);
animal = zeros(size(normptsh,1),1);
% celltypevector = zeros(size(normptsh,1),1);
for n = 1:size(normptsh,1)
    brainloc(n,1) = AreaDate{n,3};
    blackout(n,1) = str2double(AreaDate{n,8});

    if strcmp(AreaDate{n,7},'Magneto')==1
        animal(n,1) = 1;
    elseif strcmp(AreaDate{n,7},'Tacoma')==1
        animal(n,1) = 2;
    end
    
%     for m = 1:length(celltype)
%     if AreaDate{n,6}==celltype(m,1)
%         celltypevector(n,1) = celltype(m,2);
%     end
%     end
end

% Identify neurons of interest
neurons = find(normptsh(:,1)==1 & brainloc(:)==1 & blackout(:)<=3 & animal(:)>0);
% ^^^^^^^^^IMPORTANT^^^^^^^^^^^^
% For stim, normptsh(:,1)==1; For sham, normptsh(:,1)==0 

% Isolate activity profiles of interest, in this case, first 150 ms after TMS pulse
waves=normptsh(neurons,3+gauss_size+1+Tstart:3+gauss_size+1+Tend);

% K means clustering
% Evaluate range of K values to find optimal K
% Optimal is defined by the chosen criterion (CalinskiHarabasz, DaviesBouldin, silhouette, gap)
kmax = 8;
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

%% Plot results
counter=0;
idx_count = zeros(numgroups,1);
pop_psth = [];
percent_inten_cluster = zeros(numgroups,9);
percent_type_cluster = zeros(numgroups,5);

words = 'M1 Stim: All';
% interest = 1; % change this value to indicate cell type of interest!

for n=1:numgroups
pos = find(idx==n);% & celltypevector(neurons)==0); % change this value to indicate cell type of interest! 
if length(pos)<2
    continue
end
counter=counter+1;
cols = 6;

% Count n for each cluster
idx_count(n) = length(pos);

figure(3)
subplot(4,cols,counter)
hold on
imagesc(normptsh(neurons(pos),3+gauss_size:end-gauss_size),[-1 1]);
line([501 501],[0 length(pos)+1],'Color','k')
if counter==ceil(cols/2)
title(words)
end
axis([0 1000 0 length(pos)+1])

pop_psth(n,:) = nanmean(normptsh(neurons(pos),3+gauss_size:end-gauss_size));
subplot(4,cols,counter+cols)
hold on
plot(pop_psth(n,:),'r')
plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
line([tbase+1 tbase+1],[-1 1],'Color','k')
axis([0 1000 -0.5 1])
% title(['Cluster ' num2str(n)])

for j = 1:9
    intens = length(find(normptsh(neurons,2)==j*10));% & celltypevector(neurons)==0)); % change this value to indicate cell type of interest!
    inten_cluster = length(find(normptsh(neurons(pos),2)==j*10));
    percent_inten_cluster(counter,j) = 100*inten_cluster./intens;
end

subplot(4,cols,counter+2*cols)
hold on
bar(10:10:90,percent_inten_cluster(counter,:))
axis([0 100 0 100])
ylabel('Percentage') % percent of all blocks at this intensity that assign to this cluster
if counter==ceil(cols/2)
title('Percentage of These Cell Types at This Intensity that Assigns to This Cluster')
end

% for c = 1:5
%    types = length(find(celltypevector(neurons)==(c-1)));
%    type_cluster = length(find(celltypevector(neurons(pos))==(c-1)));
%    percent_type_cluster(counter,c) = 100*type_cluster./types;
% end

% subplot(4,cols,counter+3*cols)
% hold on
% bar(0:1:4,percent_type_cluster(counter,:))
% axis([-0.5 4.5 0 100])
% ylabel('Percentage') % percent of all blocks of this cell type that assign to this cluster
% if counter==ceil(cols/2)
% title('Percentage of These Cell Types that Assigns to This Cluster, 0=Unclassified, 1=Axon, 2=Inhibitory, 3=Excitatory, 4=Other')
% end
end

idx_count