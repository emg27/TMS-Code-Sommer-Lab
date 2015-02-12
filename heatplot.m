clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
close all

tbase=500; %amount of baseline collected
ta=500; %amount of time after the TMS pulse
gauss_size=30;

counter=0; %Will count the number of cells used 
allptsh=[];
for k=1:size(s,2)
    if(length(s(k).Pulses)>0) & median(diff(s(k).Pulses))>4
        pulses=s(k).Pulses;
        firerate=s(k).FireRate;
        for g=1:max(s(k).clusters)
            cluster=find(s(k).clusters==g);
            if length(cluster)>0 & length(cluster)>length(pulses)
                counter=counter+1;
                psth1block(pulses,tbase+gauss_size,ta+gauss_size, 1000*s(k).times(cluster), gauss_size,0)
            end
        end
    end
end