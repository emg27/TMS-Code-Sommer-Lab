%% Prepare Space
% clear; clf; clc;

% Notes:
% The purpose of this script is to create the PSTHs and correlation plots.
% This script only requires an s struct


%% Parameters
% % Control parameters of Interest
% % s_proxy = s(sSubInds);
% % Get exitatotry and inhibitory indices
% waveName = char(zeros([length(waveCellStruct), 8]));
% waveBlock = zeros([length(waveCellStruct), 1]);
% waveID = zeros([length(waveCellStruct), 1]);
% waveType = zeros([length(waveCellStruct), 1]);
% waveFullName = char(zeros([length(waveCellStruct), 10]));
% cellType = 5;
% cellName = '';
% 
% if(cellType == 1)
%     cellName = 'Axon';
% elseif(cellType == 2)
%     cellName = 'Inhibitory';
% elseif(cellType == 3)
%     cellName = 'Excitatory';
% elseif(cellType == 4)
%     cellName = 'Unknown';
% elseif(cellType == 5)
%     cellName = 'Unfit';
% end
% 
% for i = 1:length(waveCellStruct)
%     proxyName = strsplit(waveCellStruct(i).Name, '_');
%     waveName(i,:) = proxyName{1};
%     waveBlock(i,:) = str2num(proxyName{2});
%     waveID(i) = waveCellStruct(i).ID(end);
%     waveType(i) = waveCellStruct(i).CellType;
%     waveFullName(i,:) = [waveName(i,:) '_' num2str(waveID(i))];
% end
% 
% sID = nan([1, length(s)]);
% for i = 1:length(s)
% % for i = 50
%     proxyName = strsplit(s(i).Name, '_');s
%     proxyOne = proxyName{1};
%     proxyBlock = str2num(proxyName{2});
%     indProxyOne = strmatch(proxyOne, waveName, 'exact');
%     indProxyBlock = find(waveBlock == proxyBlock);
%     indProxy = intersect(indProxyOne, indProxyBlock);
%     if(length(indProxy) > 0)
%         sID(i) = waveType(indProxy(1));
%     end
% end


% File Parameters
s_proxy = s;
stimshamlist = zeros([length(s_proxy), 1]);
intensitylist = zeros([length(s_proxy), 1]);
stimsham = 'Sham';

% Analysis Parameters
intlist = ['10';'20';'30';'40';'50';'60';'70';'80';'90']; % Intensities to Plot
FR = 25000;             % Max FR in files, SPS (Hz)
gauswidth = 20;         % ms: Total numerical width of gaussian
gausstd = 3;            % ms: STD of Gaussian
plotoff = gauswidth*3;  % ms: #ms of data you want to count outside of plot window
gaus = normpdf([-(gauswidth/2*FR/1000):(gauswidth/2*FR/1000)], 0, gausstd*FR/1000);
gaus = gaus/max(gaus);
tmin = 200 + plotoff;
tmax = 1500 + plotoff;
matlen = round(((tmin + tmax) / 1000 * FR))+1;
zeroproxy = round((tmin/1000 * FR)+1);
timeind = [-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)];


%% Collect relevant parameters
if(strcmp(stimsham, 'Stim'))
    stimshamnum = 1;
elseif(strcmp(stimsham, 'Sham'))
    stimshamnum = 0;
elseif(strcmp(stimsham, 'Stim (COIL MOVED)'))
    stimshamnum = 2;
elseif(strcmp(stimsham, 'Sham (COIL MOVED)'))
    stimshamnum = 3;
else
    stimshamnum = -1;
end


%% Calculate
% Find indeces of desired plot parameters (i.e config, intensity)
for i = 1:length(s_proxy)
    intProx = str2num(s_proxy(i).Intensity{1});
    if(length(intProx) > 0)
        intensitylist(i) = intProx;
    end
    if(strcmp(s_proxy(i).Stim{1}, 'Stim'))
        stimshamlist(i) = 1;
    elseif(strcmp(s_proxy(i).Stim{1}, 'Sham'))
        stimshamlist(i) = 0;
    elseif(strcmp(s_proxy(i).Stim{1}, 'Stim (COIL MOVED)'))
        stimshamlist(i) = 2;
    elseif(strcmp(s_proxy(i).Stim{1}, 'Sham (COIL MOVED)'))
        stimshamlist(i) = 3;
    else
        stimshamlist(i) = -1;
    end
end

% Store average pre/post pulse PSTH of one file/block (conv of all spikes / # Pulses)
convproxyblockavg = nan([length(s_proxy), length([-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)])]);
% Store average PSTH across each intensity (sum(PSTHs at intensity i) / # files
convproxysum = zeros([9, length([-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)])]);
% Stores baseline-subtracted PSTH for every single pulse that ever existed
convproxypulse = nan([length(s_proxy)*20, length([0 : (1000/FR) : (tmax + gauswidth)])]);
% Number of blocks at each intensity
blockcount = zeros([1, size(intlist, 1)]);
% Index of files in s-struct in the order of intensitylist
proxyindall = [];
convpulseind = 1;
for k = 1:length(intlist)
    intense = intlist(k,:);

    % Get total number of pulses that meet criteria
    proxyind = find(stimshamlist == stimshamnum & intensitylist == str2num(intlist(k,:)));
    blockcount(k) = length(proxyind)
    totpulse = 0;
    for i = 1:length(proxyind)
        totpulse = totpulse + length(s_proxy(proxyind(i)).Pulses);
    end

    % Get all times
    for i = 1:length(proxyind)
        pulse = s_proxy(proxyind(i)).Pulses;
        times = s_proxy(proxyind(i)).times((s_proxy(proxyind(i)).clusters ~= 0) & (s_proxy(proxyind(i)).clusters ~= 9)) * 1000;
        proxytimes = zeros([length(pulse), matlen]);
        for j = 1:length(pulse)
            proxytimes_b = times(...
                (times > (pulse(j) - tmin)) &...
                (times < pulse(j))) - pulse(j);
            proxytimes_a = times(...
                (times > pulse(j)) &...
                (times < (pulse(j) + tmax))) - pulse(j);
            proxytimes(j,round(proxytimes_b/1000 * FR) + zeroproxy) = 1;
            proxytimes(j,round(proxytimes_a/1000 * FR) + zeroproxy) = 1;
            convproxy = conv(proxytimes(j,:), gaus);
            convproxypulse(convpulseind,:) = convproxy(zeroproxy:end) - mean(convproxy(1:zeroproxy));
            convpulseind = convpulseind + 1;
        end
        % Avg convolution per pulse
        convproxyblockavg(proxyind(i),:) = (conv(sum(proxytimes,1), gaus) / length(pulse));             
    end
    % Avg convolution per pulse per file
    convproxysum(k,:) = sum(convproxyblockavg(proxyind,:),1,'omitnan') / (sum(~isnan(proxyind)));       
    k
    proxyindall = [proxyindall; proxyind];
end

% LateInd = [33   236   549   559   562   576   585   603   646   662   672   703   708   720   721   724   784   846]
% EarlyInd = [11   194   213   226   252   301   316   331   344   354   362   370   373   377   389   413   439   453   469   474   480   487   501   535   542   553   633   716   221   446   475   730   739   749   834]


%% Calculate Baslines
animal = 'Population'
% PSTHs for grouping data into high and low intensities
desindgrouplow = find(stimshamlist == stimshamnum & ismember(intensitylist, [10 20 30 40]));
desindgrouphigh = find(stimshamlist == stimshamnum & ismember(intensitylist, [60 70 80 90]));
convproxygroup = mean(convproxyblockavg(desindgrouplow,:),1);
convproxygroup(2,:) = mean(convproxyblockavg(desindgrouphigh,:),1);
convproxygroupnorm = convproxygroup;

% PSTHs for grouping data into regular intensities (i.e. 10, 20, 30...)
convproxynorm = convproxysum;

% Subtract baseline PSTH value, no normalization
convproxybase = mean(convproxysum(:,1:zeroproxy),2,'omitnan');
convproxygroupbase = mean(convproxygroup(:,1:zeroproxy),2,'omitnan');
for i = 1:length(convproxybase)
    convproxynorm(i,:) = convproxysum(i,:) - convproxybase(i);
end
for i = 1:length(convproxygroupbase)
    convproxygroupnorm(i,:) = convproxygroup(i,:) - convproxygroupbase(i); 
end

% convproxyblockbase = mean(convproxyblockavg(proxyindall,1:zeroproxy),2);
% convproxyblocknorm = convproxyblockavg(proxyindall,:);
% for i = 1:length(convproxyblockbase)
%     convproxyblocknorm(i,:) = convproxyblockavg(proxyindall(i),:) - convproxyblockbase(i,:);
% end


%% Plot
% Figure Parameters
titleSize = 42;
labelSize = 32;
legendSize = 24;
tickSize = 32;  
% Which Figure to Plot
convplot = convproxynorm;


figure(1);
clf; 
c = colormap(jet(size(convplot,1)));
for ii = 1:(size(convplot,1))
    hold on
%     errorbar(timeind, convplot(ii,:), convproxySTE(ii,:),'Color',c(ii,:)','LineWidth',1);
      plot(timeind, convplot(ii,:),'Color',c(ii,:)','LineWidth',3);
%       plot(timeind, convplot(ii,:) + convproxySTE(ii,:),'Color',c(ii,:)','LineWidth',1);
end
set(gca,'FontSize',tickSize)
hold off
% [legh,objh,outh,outm] = legend(intlist);
% set(objh,'LineWidth',2);
AX = legend([intlist repmat(' | n=', [size(intlist,1),1]) num2str(blockcount')]);
AX.FontSize = legendSize;
xlim([-tmin+plotoff, tmax-plotoff]);
ylim([min(min(convplot)), max(max(convplot))]);
xlabel('Time (ms)', 'FontSize', labelSize);
ylabel('Average Spike Density', 'FontSize', labelSize);
title(['Response to ' stimsham ' TMS'], 'FontSize', titleSize);% Not Significant ' stimsham], 'FontSize', titleSize);
hold on;
line([0,0],[0,1], 'LineWidth', 2, 'Color', 'k');
hold off;


figure(2)
colormap(parula(100));
h = imagesc(str2num(intlist), str2num(intlist), corr(convplot'), [-.5,1]);
colorbar;
set(gca,'FontSize',tickSize)
xlabel('Intensity', 'FontSize', labelSize);
ylabel('Intensity', 'FontSize', labelSize);
title(['Correlation Coefficients Between Intensities']);
% saveas(gcf, ['CorrPlots\' alist{k} '.jpg']);



figure(3);
clf;
colormap(jet(100));
% Heatplot of same graph
intplot = [90 80 70 60 50 40 30 20 10];
intplot = intplot(end:-1:1);
% intplot = intplot(end:-1:1);
imagesc(timeind, intplot, convproxynorm);
xlabel('Time (ms)');
ylabel('Stimulation Intensity (% Machine Max)');
title([animal ' HeatMap of Spike Density vs. Intensity ' stimsham]);

