% Hierarchical clustering
% Number of clusters is determined after viewing dendrogram of clusters

% clear
% [filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
% load([pathname filename])
% close all
clf

Tstart = 500;
Tend = Tstart+150; 

brainloc = zeros(length(normptsh),1);
blackout = zeros(length(normptsh),1);
animal = zeros(length(normptsh),1);
for n = 1:length(normptsh)
    brainloc(n,1) = AreaDate{n,3};
    blackout(n,1) = str2double(AreaDate{n,8});
    if strcmp(AreaDate{n,7},'Magneto')==1
        animal(n,1) = 1;
    elseif strcmp(AreaDate{n,7},'Tacoma')==1
        animal(n,1) = -1;
    end
end

% Identify neurons of interest
neurons=find(normptsh(:,1)==1 & brainloc(:)==1 & blackout(:)<=4 & responsive(:)==1); % responsive determined by unresponsive_blocks.m
words='M1 Stim, Blackout <= 4ms';

% Calculate dissimilarity matrix, based on euclidian distances between
% waves (of activity profiles) at each time point
waves=normptsh(neurons,3+gauss_size+1+Tstart:3+gauss_size+1+Tend);
disArray = pdist(waves,'seuclidean');

% Hierarchical clustering
cutoff = 32;                        % choose value, based on dendrogram
numgroups = 8;                      % number of clusters depends on cutoff choice
Z = linkage(disArray,'complete');
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
% To compare hierarchical vs. K means, compare "groups" vs. "IDX" (from K
% means clustering code)

figure(1)
dendrogram(Z,0,'colorthreshold',cutoff)
xlabel(words)

k1=0;
avgPSTH = [];
inten = [];
inten_cluster = [];
percent_inten_cluster = [];

for k=1:max(groups)
pos=find(groups==k);
if length(pos)<2
    continue
end
k1=k1+1;

figure(2)
subplot(3,numgroups,k1)
hold on
imagesc(normptsh(neurons(pos),3+gauss_size:end-gauss_size),[-1 1]);
line([501 501],[0 length(pos)+1],'Color','k')
title(['Group ' num2str(k)])
axis([0 1000 0 length(pos)+1])
xlabel(words)

avgPSTH(k,:) = nanmean(normptsh(neurons(pos),3+gauss_size:end-gauss_size));
subplot(3,numgroups,k1+numgroups)
hold on
plot(avgPSTH(k,:),'r')
plot(1:(tbase+ta),0*(1:(tbase+ta)),'k--')
line([tbase+1 tbase+1],[-1 1],'Color','k')
axis([0 1000 -0.5 1])
title(['Group ' num2str(k)]) 
xlabel(words)

for j = 1:9
    inten(j) = length(find(normptsh(neurons,2)==j*10));
        % calculates number of blocks at each intensity
    inten_cluster(k1,j) = length(find(normptsh(neurons(pos),2)==j*10));
        % calculates number of blocks in this cluster at each intensity
    percent_inten_cluster(k1,j) = 100*inten_cluster(k1,j)./inten(j);
end

subplot(3,numgroups,k1+2*numgroups)
bar(10:10:90,percent_inten_cluster(k1,:))
title(['Group ' num2str(k)])
axis([0 100 0 100])
ylabel('Percentage') % percent of all blocks at this intensity that assign to this cluster
xlabel(words)

end
