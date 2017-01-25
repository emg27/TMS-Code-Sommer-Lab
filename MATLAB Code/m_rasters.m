%% Prepare Space
% clear; clf; close all;
% Note: This script requires function "blockraster.m", and an s struct


%% Load relevant files
% load('SfN2015_withinfo_taco_offsetcorrect.mat');
% s_taco = s;
% load('SfN2015_withinfo_mag_offsetcorrect.mat');
% s_mag = s;
% load('SfN2015_withinfo.mat');
% load('TMS_Proxy_Taco_Named.mat');
s_taco = s;
s_mag = s;

% titleSize = 48;
% labelSize = 32;
% legendSize = 32;
% tickSize = 32;
titleSize = 20;
labelSize = 20;
legendSize = 20;
tickSize = 20;


%% Plot All Data for Stim or Sham, Across all Intensities
% t_before = 200;         % ms
% t_after = 200;          % ms
% config = 'Stim';        % Config (Stim/Sham)
% s_rast = s_taco;        % Monkey
% 
% figure(10);
% clf;
% for m = 1:9
%     subplot(3,3,m);
%     intensity = num2str(m*10);
%     rastind = 1;
%     rastcolor = 1;
%     for i = 1:length(s_rast)
%         i
%         if(strcmp(s_rast(i).Stim, config))
%             if(strcmp(cell2mat(s_rast(i).Intensity), intensity))
%                 pulse = s_rast(i).Pulses;
%                 codes = (s_rast(i).clusters > 0);
%                 times = s_rast(i).times(codes) * 1000;
%                 % Grab window for each pulse
%                 for j = 1:length(pulse)
%                     proxytimes = times((times > (pulse(j)-t_before) & times < (pulse(j) + t_after)));
%                     for k = 1:length(proxytimes)
%                         if(rastcolor == 1)
%                             blockraster(0, proxytimes(k)-pulse(j), rastind, 'blue');
%                         else
%                             blockraster(0, proxytimes(k)-pulse(j), rastind, 'red');
%                         end
%                     end
%                     rastind = rastind + 1;
%                 end
%                 rastcolor = rastcolor+1;
%                 rastcolor = mod(rastcolor,2);
%                 rastind = rastind + 1
%             end
%         end
%     end
%     line([0 0], [0 rastind], 'Color', 'black', 'LineWidth', 1);
%     title([' M1 Response at ', intensity, '% ', config]);
%     set(gca,'FontSize',12)
%     xlabel('Time (ms)');
%     ylabel('Pulse #');
%     xlim([-t_before, t_after]);
%     ylim([0 rastind]);
% end


%% Plot Data for One Intensity, One Configuration, All files
% t_before = 200;         % ms
% t_after = 200;          % ms
% config = 'Sham';        % Config (Stim/Sham)
% intensity = '90';       % Intensity
% s_rast = s_taco;        % Monkey
% mName = 'Mags';         % Name of Monkey for Plot title
% 
% rastind = 1;            % Keep at 1
% rastcolor = 1;          % Keep at 1
% 
% figure(2);
% clf;
% for i = 1:length(s_rast)
%     if(strcmp(s_rast(i).Stim, config))
%         if(strcmp(cell2mat(s_rast(i).Intensity), intensity))
%             pulse = s_rast(i).Pulses;
%             codes = (s_rast(i).clusters > 0);
%             times = s_rast(i).times(codes) * 1000;
%             % Grab window for each pulse
%             for j = 1:length(pulse)
%                 proxytimes = times((times > (pulse(j)-t_before) & times < (pulse(j) + t_after)));
%                 for k = 1:length(proxytimes)
%                     if(rastcolor == 1)
%                         blockraster(0, proxytimes(k)-pulse(j), rastind, 'blue');
%                     else
%                         blockraster(0, proxytimes(k)-pulse(j), rastind, 'red');
%                     end
%                 end
%             end
%             str1 = s_mag(i).Name;
%             t = text(-t_before*.9,rastind,str1);
%             if(rastcolor == 1)
%                 t.Color = 'b';
%             else
%                 t.Color = 'r';
%             end
%             rastcolor = rastcolor+1;
%             rastcolor = mod(rastcolor,2);
%             rastind = rastind + 1
%         end
%     end
% end
% title([mName, ' ', intensity, ' ', config]);


%% Plot Individual Files
% Parameters
t_before = 300;         % ms 
t_after = 300;          % ms
s_rast = s_mag;         % Name of Monkey
config = 'Stim';        % Config: Stim//Sham//Stim (COIL MOVED)
exclude = [];           % List of files to exclude
code_plot = [1:8];          % Neuron ID to plot

% Get list of unique file names (aList)
for i = 1:length(s_rast)
    nproxy = s_rast(i).Name;
    n = strsplit(nproxy, '_');
    names{i} = n{1};
end
alist = unique(names);

% Select one file name
fname = '20151113';
% fname = '20140605';
% fname = '20120523'
% fname = '20150918';
% Extras
% Raster Parameters
rastind = 1;
rastcolor = 1;
totPlot = 0;
% figure(2);
% close;


figure(3);
clf;

% Search through all files in s_rast for unique dates matching criteria
% s_plot: A list of indexes in struct to plot
s_plot = s_rast(1);
plotind = 1;
for i = 1:length(s_rast)
    n_proxy = strsplit(s_rast(i).Name, '_');
    n = n_proxy(1);
    n = n{1};
    if(strcmp(n, fname))
        if(strcmp(s_rast(i).Stim, config))
            monkey = s_rast(i).Animal{1};
            n
            s_plot(plotind) = s_rast(i);
            plotind = plotind + 1;
        end
    end
end

% Sort s_plot in ascending order
lis = [];
for i = 1:length(s_plot)
    lisproxy = str2num(s_plot(i).Intensity{1});
    if(length(lisproxy) > 0)
        lis(i) = lisproxy;
    else
        lis(i) = 0;
    end
end
[B,I] = sort(lis);


% Plot raster of ordered s_list
% c = colormap(cool(size(I,2)));
c = colormap(cool(9));
intlist = [];
intind = 1;
hold on
for i = I
    s_plot(i).Stim
    if(~any(i==exclude))
        pulse = s_plot(i).Pulses;
        codes = (ismember(s_plot(i).clusters,code_plot));
        times = s_plot(i).times(codes) * 1000;
        if(sum(codes) > 0)
            intlist(intind,1) = str2num(s_plot(i).Intensity{1});
            rastcolor = round(intlist(intind,1)/10);
            intind = intind + 1;
            plotind = plotind + 1;
            % Grab window for each pulse
            for j = 1:length(pulse)
                proxytimes = times((times > (pulse(j)-t_before) & times < (pulse(j) + t_after)));
                for k = 1:length(proxytimes)
                    blockraster(0, proxytimes(k)-pulse(j), rastind, c(rastcolor,:)');
                end
            rastind = rastind + 1
            end
%             str1 = [s_plot(i).Name ' ' s_plot(i).Intensity{1} ' ' s_plot(i).Stim{1}];
            str1 = [s_plot(i).Intensity{1} '% '];
            t = text(-t_before*1.14,rastind - 10,str1);
%             t.Color = c(rastcolor,:)';
            t.FontSize = legendSize;
%             t = text(-t_before*.9,rastind,s_plot(i).Name);
            rastcolor = rastcolor + 1;
            rastind = rastind + 1;
            totPlot = totPlot + 1;
        end
    end
end
line([0 0], [0 rastind], 'Color', 'Black');
% line([-4 -4], [0 rastind], 'Color', 'Red');
% line([-20 -20], [0 rastind], 'Color', 'Red');
% line([4 4], [0 rastind], 'Color', 'Red');
% line([20 20], [0 rastind], 'Color', 'Red');

% line([-50 -50], [0 rastind], 'Color', 'Red');
% line([50 50], [0 rastind], 'Color', 'Red');
ylim([0 rastind]);
set(gca,'FontSize',tickSize)
% set(gca,'ytick',[])
set(gca,'ytick',0:rastind/totPlot:rastind);
set(gca,'yticklabel',[])
% legend(num2str(intlist), 0);
hold off;
% title([monkey, ' ', fname, ' ', config], 'FontSize', titleSize);
title([monkey ' ' fname ' Single Neuron Example for ' config ' Configuration']);
xlabel('Time (ms)', 'FontSize', labelSize);
ylabel('Intensity', 'FontSize', labelSize);


%% All Rasters
% % Get list of unique names (aList)
% for i = 1:length(s_mag)
%     nproxy = s_mag(i).Name;
%     ncode = num2str(s_mag(i).ID);
%     n = strsplit(nproxy, '_');
%     names{i} = [n{1} '_' ncode];
%     i
% end
% alist = unique(names);
% 
% % Parameters
% t_before = 300;         % ms
% t_after = 300;          % ms
% monkey = 'M1';
% rastind = 1;
% rastcolor = 1;
% totPlot = 0;
% c = colormap(jet(10));
% % figure(2);
% % close;
% 
% for k = 1:length(alist)
%     fnameProxy = strsplit(alist{k}, '_');
%     fname = fnameProxy{1};
%     fID = fnameProxy{2};
%     alist{k}
%     % Search through all files in s_mag for unique dates matching criteria
%     % s_plot: A list of indexes in struct to plot
%     s_plot_stim = s_mag(1);
%     s_plot_sham = s_mag(1);
%     plotindStim = 1;
%     plotindSham = 1;
%     for i = 1:length(s_mag)
%         n_proxy = strsplit(s_mag(i).Name, '_');
%         n = n_proxy(1);
%         n = n{1};
%         if(strcmp(n, fname))
%             if(strcmp(s_mag(i).Stim, 'Stim') & s_mag(i).ID == str2num(fID))
%                 monkey = s_mag(i).Animal{1};
%                 n
%                 s_plot_stim(plotindStim) = s_mag(i);
%                 plotindStim = plotindStim + 1;
%             end
%             if(strcmp(s_mag(i).Stim, 'Sham') & s_mag(i).ID == str2num(fID))
%                 monkey = s_mag(i).Animal{1};
%                 n
%                 s_plot_sham(plotindSham) = s_mag(i);
%                 plotindSham = plotindSham + 1;
%             end
%                 
%         end
%     end
% 
%     % Sort s_plot Stim in ascending order
%     lisStim = [];
%     lisSham = [];
%     for i = 1:length(s_plot_stim)
%         lisproxy = str2num(s_plot_stim(i).Intensity{1});
%         if(length(lisproxy) > 0)
%             lisStim(i) = lisproxy;
%         else
%             lisStim(i) = 0;
%         end
%     end
%     [B,IStim] = sort(lisStim);
%     % Sort s_plot Sham in ascneding order
%     for i = 1:length(s_plot_sham)
%         lisproxy = str2num(s_plot_sham(i).Intensity{1});
%         if(length(lisproxy) > 0)
%             lisSham(i) = lisproxy;
%         else
%             lisSham(i) = 0;
%         end
%     end
%     [B,ISham] = sort(lisSham);
%     
%     % Plot
%     figure(2);
%     clf;
%     subplot(1,2,1);
%     
%     % Plot raster of ordered s_list
%     intlist = [];
%     intind = 1;
%     hold on
%     
%     rastind = 1;
%     rastcolor = 1;
%     totPlot = 0;
%     plotind = 1;
%     for i = IStim
%         s_plot_stim(i).Stim
%         pulse = s_plot_stim(i).Pulses;
%         codes = s_plot_stim(i).clusters;
%         times = s_plot_stim(i).times * 1000;
%         if(sum(codes) > 0)
%             intstrproxy = str2num(s_plot_stim(i).Intensity{1});
%             if(length(intstrproxy) > 0)
%                 intlist(intind,1) = str2num(s_plot_stim(i).Intensity{1});
%                 rastcolor = round(intlist(intind,1)/10);
%                 intind = intind + 1;
%                 plotind = plotind + 1;
% 
%                 % Grab window for each pulse
%                 for j = 1:length(pulse)
%                     proxytimes = times((times > (pulse(j)-t_before) & times < (pulse(j) + t_after)));
%                     for l = 1:length(proxytimes)
%                         blockraster(0, proxytimes(l)-pulse(j), rastind, c(rastcolor,:)');
%                     end
%                     rastind = rastind + 1;
%                 end
%     %                 str1 = [s_plot(i).Name ' ' s_plot(i).Intensity{1} ' ' s_plot(i).Stim{1}];
%                 intvalproxy = s_plot_stim(i).Intensity{1};
%                 if(length(intvalproxy) > 0)
%                     str1 = [s_plot_stim(i).Intensity{1} '% '];
%                     t = text(-t_before*1.14,rastind - 10,str1);
%     %             t.Color = c(rastcolor,:)';
%                     t.FontSize = legendSize;
%                 end
%     %             t = text(-t_before*.9,rastind,s_plot(i).Name);
%                 rastcolor = rastcolor + 1;
%                 rastind = rastind + 1;
%                 totPlot = totPlot + 1;
%             end
%         end
%     end
%     line([0 0], [0 rastind], 'Color', 'Black');
% 
%     ylim([0 rastind]);
%     xlim([-t_before, t_after]);
%     set(gca,'FontSize',tickSize)
%     % set(gca,'ytick',[])
%     set(gca,'ytick',0:rastind/totPlot:rastind);
%     set(gca,'yticklabel',[])
%     % legend(num2str(intlist), 0);
%     hold off;
%     % title([monkey, ' ', fname, ' ', config], 'FontSize', titleSize);
%     title([alist{k} ' ' s_plot_stim(i).Animal{1} ' Stim']);
%     xlabel('Time (ms)', 'FontSize', labelSize);
%     % ylabel('Intensity', 'FontSize', labelSize);
%     
%     
%     subplot(1,2,2);
%     
%     % Plot raster of ordered s_list
%     intlist = [];
%     intind = 1;
%     hold on
%     
%     rastind = 1;
%     rastcolor = 1;
%     totPlot = 0;
%     plotind = 1;
%     for i = ISham
%         s_plot_sham(i).Stim
%         pulse = s_plot_sham(i).Pulses;
%         codes = s_plot_sham(i).clusters;
%         times = s_plot_sham(i).times * 1000;
%         if(sum(codes) > 0)
%             intstrproxy = str2num(s_plot_sham(i).Intensity{1});
%             if(length(intstrproxy) > 0)
%                 intlist(intind,1) = str2num(s_plot_sham(i).Intensity{1});
%                 rastcolor = round(intlist(intind,1)/10);
%                 intind = intind + 1;
%                 plotind = plotind + 1;
% 
%                 % Grab window for each pulse
%                 for j = 1:length(pulse)
%                     proxytimes = times((times > (pulse(j)-t_before) & times < (pulse(j) + t_after)));
%                     for l = 1:length(proxytimes)
%                         blockraster(0, proxytimes(l)-pulse(j), rastind, c(rastcolor,:)');
%                     end
%                 rastind = rastind + 1;
%                 end
%     %                 str1 = [s_plot(i).Name ' ' s_plot(i).Intensity{1} ' ' s_plot(i).Stim{1}];
%                 intvalproxy = s_plot_sham(i).Intensity{1};
%                 if(length(intvalproxy) > 0)
%                     str1 = [s_plot_sham(i).Intensity{1} '% '];
%                     t = text(-t_before*1.14,rastind - 10,str1);
%     %             t.Color = c(rastcolor,:)';
%                     t.FontSize = legendSize;
%                 end
%     %             t = text(-t_before*.9,rastind,s_plot(i).Name);
%                 rastcolor = rastcolor + 1;
%                 rastind = rastind + 1;
%                 totPlot = totPlot + 1;
%             end
%         end
%     end
%     line([0 0], [0 rastind], 'Color', 'Black');
% 
%     ylim([0 rastind]);
%     xlim([-t_before, t_after]);
%     set(gca,'FontSize',tickSize)
%     % set(gca,'ytick',[])
%     set(gca,'ytick',0:rastind/totPlot:rastind);
%     set(gca,'yticklabel',[])
%     % legend(num2str(intlist), 0);
%     hold off;
%     % title([monkey, ' ', fname, ' ', config], 'FontSize', titleSize);
%     title([alist{k} ' ' s_plot_sham(i).Animal{1} ' Sham']);
%     xlabel('Time (ms)', 'FontSize', labelSize);
%     % ylabel('Intensity', 'FontSize', labelSize);
%     saveas(gcf, ['Rasters\' alist{k} '.jpg']);
% end