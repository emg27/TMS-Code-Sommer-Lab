%%find the noise artifacts and then remove those clusters from the files.

clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])

close all

tb=1000;
ta=1000;
gauss_size=30;
last_bin=250;
base_save=int16([]);
wave_save=[];
z=1;
counter=0;
for k=1:size(s,2)
    count=k;
    if (length(s(k).Pulses)>0 ) & median(diff(s(k).Pulses)./1000)>4
        Pulses=s(k).Pulses;
        firerate=s(k).FireRate;
        for g=1:max(s(k).clusters)
            cluster=find(s(k).clusters==g);
            if length(cluster)>0
                
            figure(1);
            subplot(3,3,[1,4])
            Raster(Pulses,tb,ta,1000*s(k).times(cluster));
            title(sprintf('%s Cluster: %d',s(k).Name,g))
            axis([-tb ta 0 length(s(k).Pulses)])
            xlabel('ms')
            subplot(3,3,[2,5])
            meanWor=mean(s(k).waveforms(cluster,:));
            stdWor=std(s(k).waveforms(cluster,:));
            t1= 1000*(0:1:length(meanWor)-1)/s(k).FireRate;
            time=linspace(min(t1),max(t1),1000);
            if length(meanWor)<=1
                continue
            else
                counter=counter+1;    
            end
            meanW=spline(t1,meanWor,time);
            stdW=spline(t1,stdWor,time);
            plot_variance(time,meanW+stdW,meanW-stdW,'b');
            hold on;
            plot(time,meanW,'k');
            [width timept]=classify(meanW,firerate);
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
                base_save=[base_save temp'];
                wave_save{z,1}=meanW;
                wave_save{z,2}=time;
                z=z+1;
            end
            check=input('Is this file noise?')
            if check==1
                s(k).clusters(cluster)=0;
            end
            clf(1)
            end
        end
    end
end
[filename2 filepath]=uiputfile('C:\','Save Structure','.mat');
save([filepath filename2],'s')