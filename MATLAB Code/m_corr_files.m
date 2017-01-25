%% Prepare Space
% Standalone script;
% Outputs a struct containing 9x9 intensity cross correlations
% of the same size as the input struct;


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
         convproxy = conv(proxytimes(j,:), gaus);
         convproxypulse(j,:) = convproxy;
    end
    s_corr(i).Pulses = convproxypulse';
    i
end


%% One file correlations
% % Stim near-full
alist = ['20130808_1'; '20140124_2'; '20140529_1'; '20140530_1';...
    '20140605_1'; '20140609_1'; '20140609_2'; '20140625_1'; '20140627_1';...
    '20140627_2'; '20140630_1'; '20140630_2'; '20140701_1'; '20140701_2';...
    '20140702_1'; '20140703_1'; '20140703_2'; '20140708_1'; '20140709_1';...
    '20140709_2'; '20140710_1'; '20140710_2'; '20140902_3'; '20141009_1';...
    '20141021_1'; '20150109_1'; '20150123_1'; '20150123_2'; '20150130_1';...
    '20150203_1'; '20150730_1'; '20150806_4'; '20150806_5'; '20150810_2';...
    '20150814_3'; '20150901_2'; '20150903_3'; '20150904_4'; '20151103_1';...
    '20151106_1'; '20151113_2'; '20151204_1'; '20151209_1'];
% % Stim every other full
% alist = ['20140131_3'; '20140131_4'; '20140207_1'; '20140527_2';...
%     '20140603_1'; '20140612_1'; '20140613_1'; '20140624_1'; '20140701_3';...
%     '20140902_1'; '20150701_1'; '20150729_1'; '20150729_4'; '20150730_4';...
%     '20150731_1'; '20150828_1'; '20150901_1'; '20150904_2'; '20150904_3';...
%     '20150911_1'; '20150911_2'; '20150918_3'; '20151103_5'; '20151113_1'];
% Stim Both
% alist = ['20130808_1'; '20140124_2'; '20140529_1'; '20140530_1';...
%     '20140605_1'; '20140609_1'; '20140609_2'; '20140625_1'; '20140627_1';...
%     '20140627_2'; '20140630_1'; '20140630_2'; '20140701_1'; '20140701_2';...
%     '20140702_1'; '20140703_1'; '20140703_2'; '20140708_1'; '20140709_1';...
%     '20140709_2'; '20140710_1'; '20140710_2'; '20140902_3'; '20141009_1';...
%     '20141021_1'; '20150109_1'; '20150123_1'; '20150123_2'; '20150130_1';...
%     '20150203_1'; '20150730_1'; '20150806_4'; '20150806_5'; '20150810_2';...
%     '20150814_3'; '20150901_2'; '20150903_3'; '20150904_4'; '20151103_1';...
%     '20151106_1'; '20151113_2'; '20151204_1'; '20151209_1';...
%     '20140131_3'; '20140131_4'; '20140207_1'; '20140527_2';...
%     '20140603_1'; '20140612_1'; '20140613_1'; '20140624_1'; '20140701_3';...
%     '20140902_1'; '20150701_1'; '20150729_1'; '20150729_4'; '20150730_4';...
%     '20150731_1'; '20150828_1'; '20150901_1'; '20150904_2'; '20150904_3';...
%     '20150911_1'; '20150911_2'; '20150918_3'; '20151103_5'; '20151113_1'];

% % Sham near-full
% alist = ['20140124_2'; '20140529_1'; '20140530_1'; '20140605_1';...
%     '20140609_1'; '20140609_2'; '20140616_1'; '20140627_1'; '20140627_2';...
%     '20140630_1'; '20140630_2'; '20140708_1'; '20140709_1'; '20140709_2';...
%     '20140710_1'; '20140902_3'; '20141104_1'; '20141104_3';...
%     '20150123_1'; '20150123_2'; '20150123_4'; '20150130_1'; '20150203_1';...
%     '20150203_3'; '20150730_1'; '20150810_2'; '20150901_2'; '20151106_1';...
%     '20151113_2'; '20151204_1'];
% % Sham half-full
% alist = ['20140131_3'; '20140131_4'; '20140205_1'; '20140207_1';...
%     '20140218_2'; '20140224_1'; '20140224_2'; '20140226_1'; '20140226_2';...
%     '20140228_2'; '20140306_1'; '20140319_1'; '20140319_2'; '20140331_1';...
%     '20140404_2'; '20140416_1'; '20140430_1'; '20140509_1'; '20140612_1';...
%     '20140613_1'; '20141104_2'; '20150701_1'; '20150728_3'; '20150730_4';...
%     '20150806_4'; '20150806_5'; '20150814_3'; '20150828_1'; '20150901_1';...
%     '20150903_2'; '20150903_3'; '20150904_4'; '20150911_2'];
% % Sham Both
% alist = ['20140124_2'; '20140529_1'; '20140530_1'; '20140605_1';...
%     '20140609_1'; '20140609_2'; '20140616_1'; '20140627_1'; '20140627_2';...
%     '20140630_1'; '20140630_2'; '20140708_1'; '20140709_1'; '20140709_2';...
%     '20140710_1'; '20140902_3'; '20141104_1'; '20141104_3';...
%     '20150123_1'; '20150123_2'; '20150123_4'; '20150130_1'; '20150203_1';...
%     '20150203_3'; '20150730_1'; '20150810_2'; '20150901_2'; '20151106_1';...
%     '20151113_2'; '20151204_1';...
%     '20140131_3'; '20140131_4'; '20140205_1'; '20140207_1';...
%     '20140218_2'; '20140224_1'; '20140224_2'; '20140226_1'; '20140226_2';...
%     '20140228_2'; '20140306_1'; '20140319_1'; '20140319_2'; '20140331_1';...
%     '20140404_2'; '20140416_1'; '20140430_1'; '20140509_1'; '20140612_1';...
%     '20140613_1'; '20141104_2'; '20150701_1'; '20150728_3'; '20150730_4';...
%     '20150806_4'; '20150806_5'; '20150814_3'; '20150828_1'; '20150901_1';...
%     '20150903_2'; '20150903_3'; '20150904_4'; '20150911_2'];


corrSTD = zeros([1, size(alist, 1)]);
s_corrGood = struct('Val', []);
ASlopes = nan([2, 9, size(alist,1)]);
stimshamone = 1;
stimshamnum = 1;
for m = 1:size(alist,1)
% for m = 1
    m
    fname = strsplit(alist(m,:), '_');
    desnid = str2num(fname{2});
    fname = fname{1};
    % desnid = 2;
    % Nearly Full fnames
    % 20140530_1: All 9 intensities, Not very clear
    % 20140605_1: Missing 60, very good
    % 20141009_1: All 9 intensities
    % 20141021_1: Missing 10 and 90
    % 20150109_1: All 9 intensities
    % 20150123_1: All 9 intensities
    % 20150123_2: All 9 intensities
    % 20150130_1: All 9 intensities, perfect! Stim and Sham
    % 20150203_1: All 9 intensities, perfect! Stim and Sham
    % 20150814_3: All 9 intensities, perfect!
    % 20151106_1: All 9 intensities, blur
    % 20151204_1: All 9 intensities, perfect!
    % 20151209_1: Missing 10

    s_proxy = s;
    desnamelist = zeros([length(s_proxy), 1]);
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

    intlist = ['10';'20';'30';'40';'50';'60';'70';'80';'90'];
    for i = 1:length(s_proxy)
        n_proxy = strsplit(s_proxy(i).Name, '_');
        n = n_proxy(1);
        n = n{1};
        nidlist(i) = s_proxy(i).ID;
        if(strcmp(n, fname))
            desnamelist(i) = 1;
        end
    end

    figure(1);
    clf;
    colormap(jet);
    desindproxy = [];
    desindproxy = find(stimshamlist == stimshamone & desnamelist == 1 & nidlist == desnid);
    colnum = 5;
    rownum = max(ceil(length(desindproxy)/colnum),2);
    s_sub = s(desindproxy);
    sigdist = zeros([1, length(desindproxy)]);
    for i = 1:size(intlist,1)
        desind = find(stimshamlist == stimshamone & intensitylist == str2num(intlist(i,:)) & desnamelist == 1 & nidlist == desnid);
        if(length(desind) > 0)
            subplot(rownum*2,colnum,i);
            A = corr(s_corr(desind(1)).Pulses);
            A(isnan(A)) = 0 ;
            Aid = (ones(size(A)) - eye(size(A)));
            A = A.*Aid;
            sigpercent = sum(sum(A)) / (length(A)*(length(A)-1));
            sigdist(i) = sigpercent;
%             imagesc(abs(A), [0,.65]);
    %         title([s_proxy(desindone(i)).Name ' ' num2str(round(sigpercent*10)/10)]);
            namestr = strsplit(s_proxy(desind(1)).Name, '_');
    %         title(['r = ' num2str(round(sigpercent*1000)/1000); namestr(2)]);
            title([num2str(intensitylist(desind(1))) ' r=' num2str(round(sigpercent*1000)/1000)], 'FontSize', 24);
            subplot(rownum*2,colnum,i+rownum*colnum);
            hist(A(find(Aid == 1)));
            xlim([-.5, 1]);
    %         ylim([0, 300]);
        end
    end
    % subplot(4,5,i+1);
    % plot(sigdist);
    % H = colorbar;
    % H.TickLabels = [0:.1:.65];
    % H.Ticks = linspace(0, 1, 11);
    % H.Limits = [0 .65];
    % H.FontSize = legendSize;
    % xlabel('Pulse 2', 'FontSize', labelSize);
    % ylabel('Pulse 1', 'FontSize', labelSize);

    s_proxy = s_sub;
    stimshamlist = zeros([length(s_proxy), 1]);
    intensitylist = zeros([length(s_proxy), 1]);
    intval = 60;
    % Non-control
    FR = 25000;             % Max FR in files, SPS (Hz)
    gauswidth = 30;         % ms
    gausstd = 3;            % ms
    plotoff = 50;
    gaus = normpdf([-(gauswidth/2*FR/1000):(gauswidth/2*FR/1000)], 0, gausstd*FR/1000);
    gaus = gaus/max(gaus);
    tmin = -50 + plotoff;
    tmax = 150 + plotoff;
    matlen = round(((tmin + tmax) / 1000 * FR))+1;
    zeroproxy = round((tmin/1000 * FR)+1);
    timeind = [-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)];

    intlist = ['10';'20';'30';'40';'50';'60';'70';'80';'90'];
    % intlist = ['50';'60';'70';'80';'90'];

    % Get Day
    % sSubInds = [];
    % for i = 1:length(s)
    %     if(length(strfind(s(i).Name, alist{56})))
    %         sSubInds = [sSubInds, i];
    %     end
    % end


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

    convproxysum = zeros([9, length([-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)])]);
    convproxyblockavg = zeros([length(s_proxy), length([-(tmin + gauswidth/2) : (1000/FR) : (tmax + gauswidth/2)])]);
    proxypulse = zeros([1, matlen]);
    convproxypulse = zeros([length(s_proxy)*20, length([0 : (1000/FR) : (tmax + gauswidth)])]);
    blockcount = zeros([1, size(intlist, 1)]);
    proxyindall = []
    convpulseind = 1;
    for k = 1:length(intlist)
        intense = intlist(k,:);

        % Get total number of pulses that meet criteria
        proxyind = find(stimshamlist == stimshamnum & intensitylist == str2num(intlist(k,:)));
        blockcount(k) = length(proxyind);
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
            convproxyblockavg(proxyind(i),:) = (conv(sum(proxytimes,1), gaus) / length(pulse));     % Avg convolution per pulse
        end
        convproxysum(k,:) = sum(convproxyblockavg(proxyind,:),1,'omitnan') / (sum(~isnan(proxyind)));              % Avg convolution per pulse per file
        proxyindall = [proxyindall; proxyind];
    end

    % Get desired indices
    desind = find(stimshamlist == stimshamnum & intensitylist == intval);
    convproxydes = convproxyblockavg(desind,:);
    convproxypulse = convproxypulse(1:(convpulseind-1),:);

    % LateInd = [33   236   549   559   562   576   585   603   646   662   672   703   708   720   721   724   784   846]
    % EarlyInd = [11   194   213   226   252   301   316   331   344   354   362   370   373   377   389   413   439   453   469   474   480   487   501   535   542   553   633   716   221   446   475   730   739   749   834]

    % Figure Parameters
    titleSize = 42;
    labelSize = 32;
    legendSize = 24;
    tickSize = 32;  

    animal = 'Population'
    % Groupings
    desindgrouplow = find(stimshamlist == stimshamnum & ismember(intensitylist, [10 20 30 40]));
    desindgrouphigh = find(stimshamlist == stimshamnum & ismember(intensitylist, [60 70 80 90]));
    convproxygroup = mean(convproxyblockavg(desindgrouplow,:),1);
    convproxygroup(2,:) = mean(convproxyblockavg(desindgrouphigh,:),1);
    convproxygroupnorm = convproxygroup;
    convproxydesnorm = convproxydes;

    % Line graphs with nice color gradients
    convproxybase = mean(convproxysum(:,1:zeroproxy),2,'omitnan');
    convproxygroupbase = mean(convproxygroup(:,1:zeroproxy),2,'omitnan');
    convproxydesbase = mean(convproxydes(:,1:zeroproxy),2,'omitnan');
    convproxynorm = convproxysum;
    for i = 1:length(convproxybase)
            convproxynorm(i,:) = convproxysum(i,:) - convproxybase(i);
    end
    for i = 1:length(convproxygroupbase)
        convproxygroupnorm(i,:) = convproxygroup(i,:) - convproxygroupbase(i); 
    end
    for i = 1:length(convproxydesbase)
        convproxydesnorm(i,:) = convproxydes(i,:) - convproxydesbase(i); 
    end

    convproxyblockbase = mean(convproxyblockavg(proxyindall,1:zeroproxy),2);
    convproxyblocknorm = convproxyblockavg(proxyindall,:);
    for i = 1:length(convproxyblockbase)
        convproxyblocknorm(i,:) = convproxyblockavg(proxyindall(i),:) - convproxyblockbase(i,:);
    end
    convplot = convproxynorm;
    c = colormap(jet(size(convplot,1)));
    Aone = zeros([9,9]);
    ASTDone = 0;
    if(length(s_sub) > 0)
        Aone = corr(convplot');
        ADiag = Aone;
        Aidone = (ones(size(Aone)) - eye(size(Aone)));
        ADiag(1:size(Aone,1)+1:end) = nan;
        
        yMin = min(ADiag(:))
        yMax = max(ADiag(:)) - yMin;
        
        for l = 1:size(Aone,1)
            yp = ADiag(l,:)';
            xp = str2num(intlist);
            xpn = xp(~isnan(yp));
            ypn = yp(~isnan(yp));
            if(length(ypn) > 0)
%                 ASlopes(:,l,m) = glmfit(xpn, (ypn-yMin)./yMax, 'binomial', 'logit');
            end
        end
        
        figure(2)
        ASTDone = std(Aone(Aidone(:) == 1), 'omitnan');
%         imagesc(str2num(intlist), str2num(intlist), Aone, [-.5,1]);
        colorbar;
        set(gca,'FontSize',tickSize)
        xlabel('Intensity', 'FontSize', labelSize);
        ylabel('Intensity', 'FontSize', labelSize);
%         title(['Correlation Coefficients Between Intensities: ' s_sub(1).Name(1:8) '_' num2str(s_sub(1).ID)]);
        title(['STD: ' num2str(ASTDone)]);
%         saveas(gcf, ['CorrPlots\' alist(m,:) '.jpg']);
        s_corrGood(m).Val = Aone;
    end
    
    figure(3)
    Aidone = (ones(size(Aone)) - eye(size(Aone)));
%     hist(Aone(Aidone(:) == 1));
    xlim([-.4, 1]);
    corrSTD(m) = std(Aone(Aidone(:) == 1), 'omitnan');
    title(['STD: ' num2str(ASTDone) ' Histogram']);
%     saveas(gcf, ['CorrPlots\' alist(m,:) '_Hist.jpg']);
%     pause;
end


%% Distributions
numcol = round(sqrt(length(s_corrGood))*1.4);
figure(1);
gausind = [-.5: .002: 1]';
gausFits1 = zeros([length(gausind), length(s_corrGood)]);
gausFits2 = zeros([length(gausind), length(s_corrGood)]);

binWidth = .1;

AICVal = zeros([2, length(s_corrGood)]);
BICVal = zeros([2, length(s_corrGood)]);
for i = 1:length(s_corrGood)
    i
    subplot(ceil(length(s_corrGood) / numcol), numcol, i);
    Aplot = s_corrGood(i).Val;
    Aplotdiag = Aplot;
    Aidone = (ones(size(Aone)) - eye(size(Aone)));
    Aplotdiag(1:size(Aplot,1)+1:end) = nan;
    Aplotdiag = Aplotdiag(~isnan(Aplotdiag));
    
%     % Calculations
    Mu1 = .5;
    Sigma1(:,:,1) = .01;
    
    Mu2 = [.2; .8];
    Sigma2(:,:,1) = .01;
    Sigma2(:,:,2) = .01;
    
    S1 = struct('mu', Mu1, 'Sigma', Sigma1);
    S2 = struct('mu', Mu2, 'Sigma', Sigma2);
    GMModel1 = fitgmdist(Aplotdiag(:),1, 'Options',statset('MaxIter',20000,'TolFun',1e-15));
    GMModel2 = fitgmdist(Aplotdiag(:),2, 'RegularizationValue',0.0001, 'Options',statset('MaxIter',202000,'TolFun',1e-15));
    muProxy1 = GMModel1.mu;
    muProxy2 = GMModel2.mu;
    sigProxy1 = [GMModel1.Sigma(1)];
    sigProxy2 = [GMModel2.Sigma(1); GMModel2.Sigma(2)];
    AICVal(:,i) = [GMModel1.AIC; GMModel2.AIC];
    BICVal(:,i) = [GMModel1.BIC; GMModel2.BIC];
    
    histogram(Aplotdiag, 'BinWidth', binWidth);
    hold on;
    gausFits1(:,i) = GMModel1.pdf(gausind);
    gausFits2(:,i) = GMModel2.pdf(gausind);
    plot(gausind, GMModel1.pdf(gausind)*length(Aplotdiag)*binWidth, 'r-');
    plot(gausind, GMModel2.pdf(gausind)*length(Aplotdiag)*binWidth, 'b-');
%     plot(gausind, gaus2_2, 'b-');
    hold off;
    title(['MixMod: ' num2str(round(GMModel1.BIC - GMModel2.BIC))]);
    xlim([-.5, 1]);
    ylim([0 20]);
end


%% Sort
% [B,Isort] = sort(max(ASlopes(2,:,:)) - min(ASlopes(2,:,:)));x
% [B,Isort] = sort(max(abs(ASlopes(2,:,:))));
[B,Isort] = sort(BICVal(1,:) - BICVal(2,:));
% [B,Isort] = sort(corrSTD);
% Plot all Ascending Heatplots
figure(2);
clf;
for i = 1:length(Isort)
    subplot(ceil(length(Isort) / numcol), numcol, i);
    Aplot = s_corrGood(Isort(i)).Val;
%     Aplot = Aplot - Aplot.*(eye(size(Aplot)));
    imagesc(str2num(intlist), str2num(intlist), Aplot);
    title(['BIC1 - BIC2: ' num2str(round(B(i)*1000)/1000)]);
end
colormap(parula(100));
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,['\bf Inter-Intensity Correlation Coefficients Sorted by Ascending BIC1 - BIC2'],'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top', 'FontSize', 24);

% All figures for logistic fit
figure(3);
clf
c = colormap(jet(9));
for i = 1:length(Isort)
% for i = 24
    subplot(ceil(length(Isort) / numcol), numcol, i);
    AProxy = s_corrGood(Isort(i)).Val;
    AProxy(1:size(AProxy,1)+1:end) = nan;
    yMin = min(AProxy(:))
    yMax = max(AProxy(:)) - yMin;
    
    hold on;
    for j = 1:size(xp,1)
        yp = AProxy(:,j);
        plot(xp, glmval(ASlopes(:,j,Isort(i)), xp, 'logit'), 'Color', c(j,:)', 'LineWidth', 1);
    end
    
%     legend(intlist);
    for j = 1:size(xp,1)
        yp = AProxy(:,j);
        plot(xp, (yp-yMin)./yMax, 'x', 'Color', c(j,:)', 'MarkerSize', 5, 'LineWidth', 1);
    end
    title(['BIC1 - BIC2: ' num2str(round(B(i)*1000)/1000)]);
%     title(num2str(B(i)));
    hold off;
end

% Re-ordered gaus fits
figure(4);
for i = 1:length(Isort)
    i
    subplot(ceil(length(s_corrGood) / numcol), numcol, i);
    Aplot = s_corrGood(Isort(i)).Val;
    Aplotdiag = Aplot;
    Aplotdiag(1:size(Aplot,1)+1:end) = nan;
    Aplotdiag = Aplotdiag(~isnan(Aplotdiag));
    
    histogram(Aplotdiag, 'BinWidth', binWidth);
    hold on;
    plot(gausind, gausFits1(:,Isort(i)).*length(Aplotdiag)*binWidth, 'r-', 'LineWidth', 2);
    plot(gausind, gausFits2(:,Isort(i)).*length(Aplotdiag)*binWidth, 'b-', 'LineWidth', 2);
%     plot(gausind, gaus2_2, 'b-');
    hold off;
    title(['BIC1 - BIC2: ' num2str(num2str(round(B(i)*1000)/1000))]);
    xlim([-1, 1]);
    ylim([0 max([gausFits1(:,Isort(i)); gausFits2(:,Isort(i))].*length(Aplotdiag)*binWidth)*1.5]);
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,['\bf Gaussian Mixture Model for Unimodal and Bimodal Fit Sorted by Ascending BIC1 - BIC2'],'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top', 'FontSize', 24);

% Distributions of BICs
figure(4);
hold on;
plot(BICVal(1,:), BICVal(2,:), '.', 'MarkerSize', 20);
hold off;
xlabel('BIC 1');
ylabel('BIC 2');



%% Single Intensity Dose-Response
corrMat = zeros([9, 9, length(s_corrGood)]);
for i = 1:length(s_corrGood)
    matProxy = s_corrGood(i).Val;
    matProxy(1:size(matProxy,1)+1:end) = nan;
    corrMat(:,:,i) = matProxy;
end

figure(5);
clf;
hold on;
c = colormap(jet(9));
corrMean = mean(corrMat,3, 'omitnan');
corrSTD = std(corrMat, 0, 3, 'omitnan');
corrSTE = corrSTD ./ sqrt(sum(~isnan(corrMat),3));
for i = 1:size(corrMean,1)
    plotind = [1:(i-1), (i+1):9];
    errorbar(str2num(intlist(plotind,:)), corrMean(plotind,i), corrSTE(plotind,i), 'Color', c(i,:)', 'LineWidth', 5);
end
hold off;
xlim([0 100]);
ylim([0 1]);
legend(intlist);
set(gca,'FontSize',tickSize)
xlabel('Intensity');
ylabel('Correlation Coefficient');
title('Inter-Intensity Correlation Coefficient Dose-Response per Intensity');