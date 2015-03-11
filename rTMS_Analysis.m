%% Load File and get Directory 
clear; clf; close all;
[filename, pathname]=uigetfile('*.mat');
load([pathname filename]);

%% Get waveforms
% Arrays to analyze
% Overall activity divided by duration of measurement
freqTotal = [];
% Activity in a given window (window size = avgWin)
windFreq = [];
% Interval between successive spikes
ISIVal = [];

% Arrays to track time
% Time of instances of each pulse
tPulse = [];
% Time of instance of each window we're averaging across
tWaveAvg = [];
% Time of instance of each wave
tWave = [];
% Time at which each ISI measurement is taken
tISI = [];


% Parameters
% Offset for looking at an entire day's worth of data.
tOff = 0;
% Window size when looking for # of cells per Window (s)
avgWin = 1;         % s
% Index of position for windFreq
windInd = 1;
% Index to plot TMS
pulseInd = [];
% Index where each file's pulse ends
pulseFileInd = [];
% Index offset between files
indOff = 0;
% Index shift in case file doesn't start at 1
shiftInd = 4;

% TMS Conditions
cond20140205 = ['20140205', '70, 5Hz, Stim', '90, 5Hz, Stim', '90, 5Hz, Stim',...
    '90, 5Hz, Sham', '90, 5Hz, Sham'];
cond20140206 = ['20140206', '70, 5Hz, Stim', '90, 5Hz, Stim', '90, 5Hz, Stim',...
    '50, 5Hz, Stim', '10, 5Hz, Stim', '90, 5Hz, Sham', '90, 10Hz, Stim'];
condArr = {cond20140205, cond20140206};


figure(1)
for j = 1:length(s)
    i = j + shiftInd;
    if(i > length(s))
        i = i - length(s);
    end
    % Find cluster(s) we want
    w1Ind = find(s(i).clusters == 1);
    tPulse = [tPulse; (s(i).Pulses)/1000 + tOff];
    ISIVal = [ISIVal; s(i).times(w1Ind(2:end)) - s(i).times(w1Ind(1:end-1))];
    tISI = [tISI; s(i).times(w1Ind(1:end-1)) + tOff];
    tWave = [tWave; s(i).times(w1Ind) + tOff];
    pulseInd = [pulseInd; repmat(min(w1Ind + indOff), size(s(i).Pulses))];
    if(length(pulseFileInd) == 0)
        pulseFileInd = length(s(i).Pulses);
    else
        pulseFileInd = [pulseFileInd; length(s(i).Pulses) + pulseFileInd(j-1)];
    end
    
    tMax = max(s(i).times);
    tNow = 0;
    % Windowed frequency
    tLook = s(i).times(w1Ind);
    while(tNow + avgWin < tMax);
        windFreq(windInd) = length(find(tLook > tNow & tLook <= tNow + avgWin))/avgWin;
        tWaveAvg(windInd) = tNow + avgWin/2 + tOff;
        tNow = tNow + avgWin;
        windInd = windInd + 1;
    end
    
    % corresponding waveforms
    wAll = s(i).waveforms;
    w1 = wAll(w1Ind, :);
    
    subplot(4,3,j);
    plot(mean(w1)', 'k.-');
    title(s(i).Name);
    freqTotal(i) = size(w1,1)/max(s(i).times);
    tOff = tOff + max(s(i).times);
    indOff = indOff + max(length(w1Ind));
end

stimPlot = repmat(max(windFreq), size(tPulse));
stimPlot = repmat(5, size(tPulse));


% Filtering
a = 1;
ord = 3;
b = repmat(1/ord, [1, ord]);
windFilt = filter(b, a, windFreq);


%% Figures
% Frequency as a function of trial
figure(2);
plot(freqTotal);
title('Frequency per file');
xlabel('trial number');
ylabel('Frequency (activity / s)');

% Frequency 
figure(3);
subplot(2,1,1);
plot(tWaveAvg, windFreq, 'k.-');
hold on;
plot(tPulse, stimPlot, 'r.');
hold off;
title('Spike density vs. Time');
xlabel('time (s)');
ylabel('Average spike count (spikes)');
subplot(2,1,2);
plot(tWaveAvg, windFilt, 'k.-');
hold on;
plot(tPulse, stimPlot, 'r.');
hold off;

% ISI
figure(4)
plot(tISI, ISIVal);
hold on;
plot(tPulse, stimPlot, 'k.');
hold off;
title('ISI vs. time');
xlabel('time (s)');
ylabel('ISI (s)');

% Wave number vs. Time
figure(5);
plotInd = 1:length(tWave);
plot(tWave, plotInd, 'b.');
hold on;
% Plot Pulses
tPulseSize = repmat(1000, size(tPulse));
plot(tPulse, pulseInd, 'r.');
line([tPulse tPulse], get(gca, 'ylim'));
% Text for Pulses
strArr = '\leftarrow';
str1 = [];
strx = tPulse(pulseFileInd);
stry = pulseInd(pulseFileInd);
text(strx, stry, str1);
hold off;
xlabel('Time (s)');
ylabel('Spike Number');
title(['Spike Count vs. Time for: ' s(1).Name]);