%% Prepare Space
% clear; clf; close all;

% Notes
% Requires s struct


%% Parameters
% Time Windows
t_b_start = 2;      % Time in ms to start looking before pulse
t_b_end = 15;       % Time in ms to end looking before pulse
t_a_start = t_b_start;      % Time in ms to start looking after pulse
t_a_end = t_b_end;       % Time in ms to end looking after pulse

% Parameters to test
intensity = [10 20 30 40 50 60 70 80 90];     % Intensities to analyze
configuration = 1;                            % 1:Stim. 0:Sham

% Logistics
% exclude = [519];     % File excluions
iter = 1:length(s);
% iter(exclude) = []; % Ultimate indeces over which we iterate

% Sample Sizes
tot_samp_size = 0;               % Running count of total # pulses
for i = 1:length(iter)
    tot_samp_size = tot_samp_size + length(s(iter(i)).Pulses);
end


%% Population-wide Counts
pop_count_a = zeros([1,tot_samp_size]);     % All spike counts after pulse
pop_count_b = zeros([1,tot_samp_size]);     % All spike counts before pulse
intens = zeros([1,tot_samp_size]);          % List of pulse intensities
config = zeros([1,tot_samp_size]);          % 1:Stim. 0:Sham
pop_ind = 1;

disp('Collecting times before and after pulse');
for i = 1:length(iter)
    pulse = s(iter(i)).Pulses;
    times = s(iter(i)).times(s(iter(i)).clusters ~= 0) * 1000;
    for j = 1:length(pulse)
        proxytimes_b = times(...
            (times > (pulse(j) - t_b_end)) &...
            (times < (pulse(j) - t_b_start)));
        proxytimes_a = times(...
            (times > (pulse(j) + t_a_start)) &...
            (times < (pulse(j) + t_a_end)));
        pop_count_b(pop_ind) = length(proxytimes_b);
        pop_count_a(pop_ind) = length(proxytimes_a);
        proxystr = str2num(s(iter(i)).Intensity{1});
        if(length(proxystr) > 0)
            intens(pop_ind) = str2num(s(iter(i)).Intensity{1});
        end
        config(pop_ind) = strcmp(s(iter(i)).Stim{1}, 'Stim');
        
        pop_ind = pop_ind + 1;
    end
end
pop_count_diff = pop_count_a - pop_count_b; % After - Before TMS
samp_count = zeros([1,length(intensity)]);  % # pulses for each intensity
for i = 1:length(intensity)
    samp_count(i) = length(find(intensity(i) == intens));
end


%% Empirical
mean_diff = zeros([2, length(intensity)]);
STE_diff = zeros([2, length(intensity)]);
% Get max length
maxstim = 0;
maxsham = 0;
for i = 1:length(intensity)
    desindstim = find(intens == intensity(i) & config == 1);
    if(length(desindstim) > maxstim)
        maxstim = length(desindstim);
    desindsham = find(intens == intensity(i) & config == 0);
    end
    if(length(desindsham) > maxsham)
        maxsham = length(desindsham);
    end
end
stimcounts = NaN([maxstim, length(intensity)]);
shamcounts = NaN([maxsham, length(intensity)]);

% Get means and STDs
for i = 1:length(intensity)
    % Stim
    desindstim = find(intens == intensity(i) & config == 1);
    mean_diff(1,i) = mean(pop_count_a(desindstim) - pop_count_b(desindstim));
    STE_diff(1,i) = std(pop_count_a(desindstim) - pop_count_b(desindstim)) / sqrt(length(desindstim));
    stimcounts(1:length(desindstim),i) = (pop_count_a(desindstim) - pop_count_b(desindstim)) *1000 / (t_a_end - t_a_start);
    % Sham
    desindsham = find(intens == intensity(i) & config == 0);
    mean_diff(2,i) = mean(pop_count_a(desindsham) - pop_count_b(desindsham));
    STE_diff(2,i) = std(pop_count_a(desindsham) - pop_count_b(desindsham)) / sqrt(length(desindsham));
    shamcounts(1:length(desindsham),i) = pop_count_a(desindsham) - pop_count_b(desindsham);
end

% g1 = num2str(intensity');
% [P, ANOVATAB, STATS] = anova1(shamcounts);
% [c,m,h,nms] = multcompare(STATS);
% ax = gca;
% ax.YTickLabel = g1(end:-1:1,:);
% xlabel('Firing Rate (Spikes/s)');
% ylabel('Intensity (% Machine Max)');
% title('Results from Tukey Post-Hoc Multiple Comparisons Test');

%% Resample: Population
disp('Resampling data randomly');
num_resamp = 10000;              % Number of resamples
pop_dist_mean = zeros([length(intensity), num_resamp]); % Resampled dist
for i = 1:length(samp_count)
    samp_size = samp_count(i);          % Sample size N at one intensity
    for j = 1:num_resamp
        samp_ind = randi([1,tot_samp_size], 1, samp_size); % N trial sample
        pop_dist_mean(i,j) = mean(pop_count_diff(samp_ind)); % Mean resamp
    end
end


%% Resample each intensity
disp('Resampling data by intensity');
samp_dist_mean = zeros([length(intensity), num_resamp]); % Resampled dist
for i = 1:length(intensity)
    inten = intensity(i);               % One intensity
    int_samp = pop_count_diff(find(...
        intens == inten &...
        config == configuration)); % Sample @ Intensity
    samp_size = samp_count(i);          % Sample size N at one intensity
    samp_int = length(int_samp);
    for j = 1:num_resamp
        samp_ind = randi([1,samp_int], 1, samp_size); % N trial sample
        samp_dist_mean(i,j) = mean(int_samp(samp_ind)); % Mean resamp
    end
end


%% Plot
clf; close all;

% Figure Parameters
titleSize = 32;
labelSize = 32;
legendSize = 32;
tickSize = 32;

% String Formatting
str_int = num2str(intensity');
str_samp = num2str(samp_count');
str_leg = [str_int, repmat(', n = ',length(str_samp),1), str_samp];
plotmin = -.2 * 1000 / (t_a_end - t_a_start);
plotmax = 1.0 * 1000 / (t_a_end - t_a_start);
titles = '';
anim = 'M1';
if(configuration == 1)
    titles = [titles 'Stim '];
elseif(configuration == 0)
    titles = [titles 'Sham '];
end
titles = [titles num2str(t_b_start) '-' num2str(t_b_end) 'ms before to '...
    num2str(t_a_start) '-' num2str(t_a_end) 'ms after '];

% Spike Rate
samp_dist_rate = samp_dist_mean * 1000 ./ (t_a_end - t_a_start); % spikes/s
pop_dist_rate = pop_dist_mean * 1000 ./ (t_a_end - t_a_start); % spikes/s
samp_dist_plot = samp_dist_rate;
pop_dist_plot = pop_dist_rate;


% Histogram distribution of mean sample differences
figure(1);
subplot(2,1,1);
hist(samp_dist_plot', 200);
xlabel('Mean Spike Rate Difference (spikes/s)');
ylabel('Frequency');
title([anim ' ' titles 'Bootstrap of Resampled Mean Sample Differences for Each Intensity']);
legend(str_leg);
xlim([plotmin, plotmax]);

% Histogram distribution of mean population differences
subplot(2,1,2);
hist(pop_dist_plot', 100);
xlabel('Mean Spike Rate Difference (spikes/s)');
ylabel('Frequency');
title('Bootstrap of Resampled Mean Population Differences for Each Intensity');
legend(str_leg);
xlim([plotmin, plotmax]);


% Bootstrapped Results, only stim
samp_se = std(samp_dist_plot');
samp_mean = mean(samp_dist_plot');
figure(2);
errorbar(intensity, samp_mean, 2*samp_se, 'LineWidth', 4, 'Color', 'Red');
hold on;
errorbar(intensity, ShamMean, 2*ShamSD, 'LineWidth', 4, 'Color', 'Blue');
hold off;
AX = legend('Stim +/- 2SE', 'Sham +/- 2SE');
AX.FontSize = legendSize;
set(gca,'FontSize',tickSize)
xlabel('Intensity (% Machine Max)', 'FontSize', labelSize);
ylabel('Firing Rate (Spikes/s)', 'FontSize', labelSize);
title([anim ' Bootstrapped Spike Firing Rate Change for ' titles], 'FontSize', titleSize);


% Empirical Results (Not Bootstrapped)
figure(2);
tdur = (t_a_end - t_a_start) / 1000;
errorbar(intensity, mean_diff(1,:)/tdur, 2*STE_diff(1,:)/tdur, 'LineWidth', 4, 'Color', 'Red');
hold on;
errorbar(intensity, mean_diff(2,:)/tdur, 2*STE_diff(2,:)/tdur, 'LineWidth', 4, 'Color', 'Blue');
hold off;
AX = legend('Stim +/- 2SE', 'Sham +/- 2SE');
AX.FontSize = legendSize;
% ylim([-20 80]);
set(gca,'FontSize',tickSize)
xlabel('Intensity (% Machine Max)', 'FontSize', labelSize);
ylabel('Firing Rate (Spikes/s)', 'FontSize', labelSize);
title([anim ' Spike Firing Rate Change from ' num2str(t_b_start)...
    '-' num2str(t_b_end) 'ms before TMS to ' num2str(t_a_start) '-' num2str(t_a_end)...
    'ms after TMS '], 'FontSize', titleSize);