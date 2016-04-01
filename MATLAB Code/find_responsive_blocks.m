% For each block, find the average +/- 1.96 SEM of baseline (pre-TMS activity).
% Then, looking at the average sdf of post-TMS activity, see if the sdf ever goes above or below those boundaries
% (i.e. >1.96*SEM above or below average baseline activity).
% If it does, it's considered TMS-responsive and is included in the clustering analysis.
% If it doesn't, it goes into "block 0, unresponsive" and we don't include it in the clustering analysis.

% Use 1.96 because it corresponds to 0.05 chance level.

%function[responsive_blocks] = find_responsive_blocks(psth_data,gauss_size)

psth_data = normptsh; % could use allptsh or normptsh

preTMS = psth_data(:,3+gauss_size+1:3+gauss_size+500);
t_end = 150; % ms after TMS pulse
postTMS = psth_data(:,3+gauss_size+500+1:3+gauss_size+500+t_end);

mean_preTMS = zeros(size(psth_data,1),1);
std_preTMS = zeros(size(psth_data,1),1);
SEM_preTMS = zeros(size(psth_data,1),1);
mean_postTMS = zeros(size(psth_data,1),1);
responsive_blocks = zeros(size(psth_data,1),1);
ind_responsive = [];
ind_unresponsive = [];

for n = 1:size(psth_data,1)
    mean_preTMS(n) = sum(preTMS(n,:))./size(preTMS,2); % calculates mean pre-TMS activity for each neuron
    std_preTMS(n) = sqrt((sum(preTMS(n,:)-mean_preTMS(n)).^2)./size(preTMS,2));
    SEM_preTMS(n) = std_preTMS(n)./sqrt(size(preTMS,2)); % calculates SEM of pre-TMS activity for each neuron
    mean_postTMS(n) = sum(preTMS(n,:))./size(postTMS,2); % calculates mean post-TMS activity for each neuron
    
    if mean_postTMS(n) > mean_preTMS(n)+3*SEM_preTMS(n) ||  mean_postTMS(n) < mean_preTMS(n)-3*SEM_preTMS(n)
    % if max(postTMS(n,:))>mean_preTMS(n)+3*SEM_preTMS(n) %||  min(postTMS(n,:))<mean_preTMS(n)-1.96*SEM_preTMS(n)
        ind_responsive = [ind_responsive,n];
        responsive_blocks(n)=1;
    else
        ind_unresponsive = [ind_unresponsive,n];
    end
end

% figure
% subplot(1,2,1); hist(psth_data(ind_unresponsive,2),9)
% title('Unresponsive')
% axis([0 100 0 400])
% subplot(1,2,2); hist(psth_data(ind_responsive,2),9)
% title('Responsive')
% axis([0 100 0 400])

for j = 1:9
    inten(j) = length(find(psth_data(:,2)==j*10));
        % calculates number of blocks at each intensity
    inten_responsive(j) = length(find(psth_data(ind_responsive,2)==j*10));
        % calculates number of responsive blocks at each intensity
    percent_responsive(j) = 100*inten_responsive(j)./inten(j);
        % calculates percent of responsive blocks at each intensity
end

figure
bar(10:10:90,percent_responsive(:))
ylabel('Percentage')

[mean_preTMS std_preTMS SEM_preTMS mean_postTMS]