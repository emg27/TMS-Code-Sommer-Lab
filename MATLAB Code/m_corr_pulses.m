%% Parameters
stimsham = 'Stim';              % Stim configuration
intval = 70;                    % Intensity to observe
% List of intensities to examine
intlist = ['10';'20';'30';'40';'50';'60';'70';'80';'90'];
sigthresh = .04;                 % Threshold of percent to count significant
FR = 25000;             % Sampling Rate in files, SPS (Hz)
% Gaussian
gauswidth = 30;         % ms
gausstd = 3;            % ms
plotoff = 50;           % Extra space to either side for plot regularity
gaus = normpdf([-(gauswidth/2*FR/1000):(gauswidth/2*FR/1000)], 0, gausstd*FR/1000);
gaus = gaus/max(gaus);
% Plot window
tmin = -50 + plotoff;
tbet = 0;
tmax = 150 + plotoff;
matlen = round(((tmin + tmax) / 1000 * FR))+1;      % Length of time window
zeroproxy = round((tmin/1000 * FR)+1);              % Index of t = 0
% Length of convolved window
timeind = [-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)];


%% Parameters to Calculate
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

s_proxy = s;
stimshamlist = zeros([length(s_proxy), 1]);
intensitylist = zeros([length(s_proxy), 1]);
nidlist = zeros([length(s_proxy), 1]);

% Find indeces of parameters
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


%% Calculate Convolution (SDF)
s_corr = struct('Pulses', []);

for i = 1:length(s_proxy)
    pulse = s_proxy(i).Pulses;
    times = s_proxy(i).times((s_proxy(i).clusters ~= 0) & (s_proxy(i).clusters ~= 9)) * 1000;
    proxytimes = zeros([length(pulse), matlen]);
    convproxypulse = zeros([length(pulse), length(timeind)]);
    
    % Look at time window for each pulse
    for j = 1:length(pulse)
         proxytimes_b = times(...
                (times > (pulse(j) - tmin)) &...
                (times < pulse(j))) - pulse(j);
         proxytimes_a = times(...
                (times > pulse(j)) &...
                (times < (pulse(j) + tmax))) - pulse(j);
         proxytimes(j,round(proxytimes_b/1000 * FR) + zeroproxy) = 1;
         proxytimes(j,round(proxytimes_a/1000 * FR) + zeroproxy) = 1;
         convproxypulse(j,:) = conv(proxytimes(j,:), gaus);
    end
    s_corr(i).Pulses = convproxypulse';
    i
end


%% Plot
% Figure Parameters
titleSize = 32;
labelSize = 32;
legendSize = 32;
tickSize = 32;

% Autocorrelate everyone and get histograms
sigpercentall = zeros([length(s_corr), 1]);
sigthreshall = zeros([length(s_corr), 1]);
s_sigdiag = struct('Ranges', []);
updown = 2;
for i = 1:length(s_corr)
    if(size(s_corr(i).Pulses,2) ~= 0)
        A = corr(s_corr(i).Pulses);
        A(isnan(A)) = 0 ;
        % Remove Diagonal
        A = A.*(ones(size(A)) - eye(size(A)));
        proxyupdown = zeros([2*updown + 1, size(A,1) - 2*updown]);
        for j = (updown+1):(size(A,1) - updown)
            proxyupdown(:,j-updown) = A((j-updown):(j+updown),j);
        end
        s_sigdiag(i).Ranges = proxyupdown;
        sigpercentall(i) = sum(sum(A)) / (length(A)*(length(A)-1));
        sigthreshall(i) = (sigpercentall(i)) > sigthresh;
    end
end

% Get max length
maxstim = 0;
maxsham = 0;
intensity = str2num(intlist);
for i = 1:length(intensity)
    desindstim = find(intensitylist == intensity(i) & stimshamlist == 1);
    if(length(desindstim) > maxstim)
        maxstim = length(desindstim);
    desindsham = find(intensitylist == intensity(i) & stimshamlist == 0);
    end
    if(length(desindsham) > maxsham)
        maxsham = length(desindsham);
    end
end
stimcounts = NaN([maxstim, length(intensity)]);
shamcounts = NaN([maxsham, length(intensity)]);
sigintmean = zeros([2, size(intlist,1)]);
sigintSTE = zeros([2, size(intlist,1)]);

for i = 1:size(intlist,1)
    % Sham
    desindsham = find(stimshamlist == 0 & intensitylist == str2num(intlist(i,:)));
    sigintmean(1,i) = mean(sigpercentall(desindsham));
    sigintSTE(1,i) = std(sigpercentall(desindsham))/sqrt(length(desindsham));
    shamcounts(1:length(desindsham),i) = sigpercentall(desindsham);
    % Stim
    desindstim = find(stimshamlist == 1 & intensitylist == str2num(intlist(i,:)));
    sigintmean(2,i) = mean(sigpercentall(desindstim));
    sigintSTE(2,i) = std(sigpercentall(desindstim))/sqrt(length(desindstim));
    stimcounts(1:length(desindstim),i) = sigpercentall(desindstim);
end


%%
% Figure Parameters
titleSize = 48;
labelSize = 32;
legendSize = 24;
tickSize = 32;

desindone = find((stimshamlist == stimshamnum) & (intensitylist == intval));

subcol = round(sqrt(length(desindone))*1.4);           % Num cols in subplot
subrow = ceil(length(desindone)/ subcol);              % Num rows in subplot

% Percentage significant vs. intensity
sigfinalstim = zeros([length(intlist), 1]);
sigfinalsham = zeros([length(intlist), 1]);
for i = 1:length(intlist)
    desindstim = find(intensitylist == str2num(intlist(i,:)) & stimshamlist == 1);
    desindsham = find(intensitylist == str2num(intlist(i,:)) & stimshamlist == 0);
    sigfinalstim(i) = sum(sigthreshall(desindstim)) / length(desindstim) * 100;
    sigfinalsham(i) = sum(sigthreshall(desindsham)) / length(desindsham) * 100;
end

% One intensity all correlations
figure(1);
clf;
colormap(parula);
for i = 1:length(desindone)
    subplot(subrow, subcol, i);
    A = corr(s_corr(desindone(i)).Pulses);
    A(isnan(A)) = 0 ;
    A = A.*(ones(size(A)) - eye(size(A)));
    imagesc(abs(A), [0,1]);
end