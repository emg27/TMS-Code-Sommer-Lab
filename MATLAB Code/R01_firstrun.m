%% Load File and get Directory 
% filename='Jessica_data.mat';
clear
% .mat file contains struct "s" with "n" spikes in the following parameters:
%   FireRate
%   Pulses (TMS Pulses)
%   times (time in seconds of each recorded spike)
%   clusters (cluster that each of the n spikes belongs to)
%   waveforms (waveforms of each of the n spikes)
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
% filepath=uigetdir('C:\','Pick Where is save the images');
% close all
%s=rmfield(s,'width');
%s=rmfield(s,'timept');

tb=100; % time before, 100ms before
ta=100; % time after, 1000ms after
gauss_size=30;
last_bin=250;
base_save=int16([]);
wave_save=[];
z=1;
counter=0;
alignment=[]; %Check if the TMS Pulses are aligned correctly
for k=1:size(s,2)
    count=k;
    if (length(s(k).Pulses)>0 ) & median(diff(s(k).Pulses)./1000)>4%&  size(s(k).Intensity,1)>0 &... length(s(k).Pulses)<=25 &
            %strcmp(s(k).Stim(1),'Stim')==1) %& mean(diff(s(k).Pulses)./1000)<4)% & mean(diff(s(k).Pulses)./1000)>4)
        Pulses=s(k).Pulses;%+16.92*ones(size(s(k).Pulses));
        firerate=s(k).FireRate;
        % Iterate through clusters and generate figures
        for g=1:max(s(k).clusters)
            cluster=find(s(k).clusters==g);
            if length(cluster)>0 %& length(cluster)~=length(s(k).Pulses)
            figure%(1);
            
            subplot(3,3,[1,4])
            Raster(Pulses,tb,ta,1000*s(k).times(cluster));
            title(sprintf('%s Cluster: %d',s(k).Name,g))
            axis([-tb ta 0 length(s(k).Pulses)])
            xlabel('ms')
            subplot(3,3,[2,5])
            meanWor=mean(s(k).waveforms(cluster,:));
            stdWor=std(s(k).waveforms(cluster,:));
            t1= 1000*(0:1:length(meanWor)-1)/s(k).FireRate; %linspace(-1,1.5,length(meanW));%
            time=linspace(min(t1),max(t1),1000);
            if length(meanWor)<=1 %| t1(end)<2
%                 if t1(end)<2
%                     close
%                 end
                continue
            else
                counter=counter+1;    
            end
            meanW=spline(t1,meanWor,time);
            stdW=spline(t1,stdWor,time);
            plot_variance(time,meanW+stdW,meanW-stdW,'b');
            hold on;
            plot(time,meanW,'k');
            %check=input('Spike?:\n')
            %if check==1
            [width timept]=classify(meanW,time);%firerate);
%             [s(k).width s(k).timept]=classify(meanW,firerate);
%                  if mean(isnan(s(k).timept))<1
%                      plot(time(s(k).timept),meanW(s(k).timept),'go')
%                  end
            s(k).width(g)={width};
            s(k).timept(g)={timept};
            if mean(isnan(timept))<1
            plot(time(timept),meanW(timept),'go')
            end
            title(sprintf('Average Waveform- Total Number of Spikes: %d',length(s(k).times(cluster))))
            xlabel('ms')
            xlim([min(time) max(time)])
            subplot(3,3,[7:9])
            blockraster(Pulses,1000*s(k).times(cluster),0,[0 0 1])
            xlabel('Entire Block Raster')
            axis([0 1000*s(k).times(end) -2 2])
            subplot(3,3,[3,6])
            isi=isigraph(1000*s(k).times(cluster),0,1000*s(k).times(cluster(end)),1,last_bin);
            title(['Position in Structure: ' num2str(k)])
            xlabel([num2str(100*sum(isi(1:3))/sum(isi)) '% multiunit activity in first 3ms ISI bins'])
            xlim([0 100])
            subplot(3,3,[7:9])
            xlabel([num2str(100*sum(isi(1:2))/sum(isi)) '% multiunit activity in first 2ms ISI bins'])
            if 100*isi(1)/sum(isi)<8 & ~isnan(100*isi(1)/sum(isi)) & size(isi(isi~=0),1)>1
                temp(1)=k;
                temp(2)=g;
               % temp(3)=str2num(s(k).Intensity{1});
                base_save=[base_save temp'];
                wave_save{z,1}=meanW;
                wave_save{z,2}=time;
                z=z+1;
            end
            if g==9 | length(cluster)==length(s(k).Pulses)
                if length(s(k).Pulses)==length(s(k).times(cluster))
                    diffPulse=s(k).Pulses-1000*s(k).times(cluster);
                elseif length(s(k).Pulses)<length(s(k).times(cluster))
                    diffPulse=s(k).Pulses-1000*s(k).times(cluster(1:length(s(k).Pulses)));
                elseif length(s(k).Pulses)>length(s(k).times(cluster))
                    diffPulse=s(k).Pulses(1:length(cluster))-1000*s(k).times(cluster);
                else
                    continue
                end
                alignment(k,:)=[k mean(diffPulse)];
            end
            close
            %plot(time,s(k).waveforms(cluster,:))
            %title('All Waveforms')
%             figure
%             psth1block(s(k).Pulses,tb+gauss_size,ta+gauss_size,1000*s(k).times(cluster),gauss_size,0);
%             title(sprintf('%s Cluster: %d',s(k).Name,g))
%             xlim([-tb ta])
% if counter==1
%     pause
% end
% print(gcf,'-dpng', [filepath '\spike' num2str(counter)])          
% clf(1)
            end
        end
    else
%    index=regexp(s(k).Name,'_B');
%     if size(index)<=0
%         index=regexp(s(k).Name,'_C');
%     end
%     if size(index)>0
%         for p=1:max(s(k).clusters)
%             code=find(s(k).clusters==p);
%             if length(code)>0
%                 figure
%                 subplot(3,3,[1,4])
%                 Raster(1000*linspace(1,s(k).times(end),15),100,100,1000*s(k).times(code))
%                 title(sprintf('%s Cluster: %d',s(k).Name,p))
%                 axis([-100 100 0 15])
%                 subplot(3,3,[2,5])
%                 meanW=mean(s(k).waveforms(code,:));
%                 stdW=std(s(k).waveforms(code,:));
%                 time=(0:1:length(meanW)-1)/50;
%                 plot_variance(time,meanW+stdW,meanW-stdW,'b')
%                 hold on;
%                 plot(time,meanW,'k');
%                 [s(k).width s(k).timept]=classify(meanW,firerate);
%                 if mean(isnan(s(k).timept))<1
%                     plot(time(s(k).timept),meanW(s(k).timept),'go')
%                 end
%                 title(sprintf('Average Waveform- Total Number of Spikes: %d',length(s(k).times(code))))
%                 subplot(3,3,[7:9])
%                 blockraster(linspace(0,s(k).times(end),15)',s(k).times(code),0,[0 0 1])
%                 xlabel('Entire Block Raster')
%                 subplot(3,3,[3,6])
%                 plot(time,s(k).waveforms(code,:))
%                 title('All Waveforms')
%             end
%         end
%     end
    end
end

%save([filename],'s')