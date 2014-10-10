% filename='Jessica_data.mat';
% load(filename)
%close all
%s=rmfield(s,'width');
%s=rmfield(s,'timept');
tb=100;
ta=500;
gauss_size=30;
last_bin=250;
for k=2%1:size(s,2)
    if length(s(k).Pulses)>0
        s(k).Pulses=s(k).Pulses-16.92*ones(size(s(k).Pulses));
        firerate=s(k).FireRate;
        for g=1:max(s(k).clusters)
            cluster=find(s(k).clusters==g);
            if length(cluster)>0
            figure;
            subplot(3,3,[1,4])
            Raster(s(k).Pulses,tb,ta,1000*s(k).times(cluster));
            title(sprintf('%s Cluster: %d',s(k).Name,g))
            axis([-tb ta 0 length(s(k).Pulses)])
            subplot(3,3,[2,5])
            meanW=mean(s(k).waveforms(cluster,:));
            stdW=std(s(k).waveforms(cluster,:));
            t1= linspace(-1,1.5,length(meanW));%1000*(0:1:length(meanW)-1)/s(k).FireRate;
            time=linspace(-1,1.5,1000);
            meanW=spline(t1,meanW,time);
            stdW=spline(t1,stdW,time);
            plot_variance(time,meanW+stdW,meanW-stdW,'b');
            hold on;
            plot(time,meanW,'k');
            %check=input('Spike?:\n')
            %if check==1
            [width timept]=classify(meanW,firerate);
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
            xlim([-.8 1.3])
            subplot(3,3,[7:9])
            blockraster(s(k).Pulses,1000*s(k).times(cluster),0,[0 0 1])
            xlabel('Entire Block Raster')
            axis([0 1000*s(k).times(end) -2 2])
            subplot(3,3,[3,6])
            isigraph(1000*s(k).times(cluster),0,1000*s(k).times(cluster),5,last_bin);
            xlim([0 last_bin])
            %plot(time,s(k).waveforms(cluster,:))
            %title('All Waveforms')
            figure
            psth1block(s(k).Pulses,tb+gauss_size,ta+gauss_size,1000*s(k).times(cluster),gauss_size,0);
            title(sprintf('%s Cluster: %d',s(k).Name,g))
            xlim([-tb ta])
            end
        end
    else
    index=regexp(s(k).Name,'_B');
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