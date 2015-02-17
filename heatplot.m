clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
close all

tbase=500; %amount of baseline collected
ta=500; %amount of time after the TMS pulse
gauss_size=5;

counter=0; %Will count the number of cells used 
allptsh=[];
normptsh=[];
for k=1:size(s,2)
    if(length(s(k).Pulses)>0) & median(diff(s(k).Pulses))>4 & ...
            size(s(k).Stim,1)>0 & strcmp(s(k).Stim(1),'Stim')==1 &...
            size(s(k).Intensity,1)>0 & strcmp(s(k).Intensity(1),'90')==1
        pulses=s(k).Pulses;
        firerate=s(k).FireRate;
        for g=1:max(s(k).clusters)
            cluster=find(s(k).clusters==g);
            if length(cluster)>0 & length(cluster)>length(pulses)
                counter=counter+1;
                [spk_d,trl_fr,bin_start_times,baseline,mean_trl_fr]=...
                    psth1block(pulses,tbase+gauss_size,ta+gauss_size, 1000*s(k).times(cluster), gauss_size,0);
                close;
                allptsh=[allptsh; mean_trl_fr];
                normptsh=[allptsh; mean_trl_fr/max(abs(mean_trl_fr))];
            end
        end
    end
end
figure
%subplot(1,2,1)
imagesc(allptsh)
line([500 500],[0 828],'Color','k')
ylabel('Cells-All Intensities')
xlabel('Time (ms)')
title('Stim')
% figure
% %ubplot(1,2,2)
% imagesc(normptsh,[-.01,.01])
% line([500 500],[0 828],'Color','k')
% ylabel('Cells-All Intensities')
% xlabel('Time (ms)')
%title('Norm Sham')