%%analysisR01: Initial Run Through of analysis code basically uploads all
%%the data from the Spike2 .mat files into a single .mat file which is make
%%it easier to access the block information. 
clear; close all;

%fileloc: This is the location of the .mat files, change this folder
%location to where the .mat files are.
fileloc=['Jessica\'] %Change to where you stored the .mat files
newfilename='Jessica_data'; %Change to the desired name you want save with.
saveloc=['']; %Change to the location where you want to save to
listing=dir([fileloc '*.mat']); 

if exist([newfilename '.mat'],'file')
    load([newfilename '.mat']);
    g=size(s,1);
else
    s=struct('Name',{},'FireRate',[],'Pulses',{},'times',{},'clusters',{},'waveforms',{},'Stim',{},'Intensity',{},'Area',{},'Depth',{})
    g=0;
end
for k=g+1:g+size(listing,1);
    file=listing(k).name(1:end-4);
    load([fileloc file '.mat'],'-regexp', '_Ch');
    %I'm going to assume that all waveforms are in Ch2 or Ch3
    if exist(['V' file '_Ch4'])
       eval(['data=V' file '_Ch4'])
    %%Save data into the structure
    s(k).Name=file;
    s(k).FireRate=1/data.interval;
    s(k).times=data.times;
    s(k).clusters=data.codes(:,1);
    s(k).waveforms=data.values;
    elseif exist(['V' file '_Ch3'])
        eval(['data=V' file '_Ch3'])
            %%Save data into the structure
    s(k).Name=file;
    s(k).FireRate=1/data.interval;
    s(k).times=data.times;
    s(k).clusters=data.codes(:,1);
    s(k).waveforms=data.values;
    end
    if exist(['V' file '_Ch32'])
        eval(['time_trig=V' file '_Ch32.times']);
        time_trig=time_trig+0.0197*ones(length(time_trig),1);
        s(k).Pulses=1000*time_trig;
    else
        figure(1)
        for b=max(data.codes(:,1)):-1:1
            subplot(2,5,double(b))
            plot(data.values(data.codes(:,1)==b,:)')
            title(sprintf('Code is %d',b))
        end
        for m=max(data.codes(:,1)):-1:1
            figure(2)
            plot(data.values(data.codes(:,1)==m,:)')
            check=input('Is this the TMS pulse?\n');
            if check==1
                s(k).Pulses=1000*data.times(data.codes(:,1)==m)
                sprintf('Exit \n')
                clf(1)
                break
            end
            sprintf('Still Here \n')
        end
        clf
    end
    clear -regexp '_Ch'
end

save([saveloc newfilename '.mat'],'s')