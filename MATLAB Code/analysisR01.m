%%analysisR01: Initial Run Through of analysis code basically uploads all
%%the data from the Spike2 .mat files into a single .mat file which is make
%%it easier to access the block information.
clear; close all;

%fileloc: This is the location of the .mat files, change this folder
%location to where the .mat files are.
[fileloc]=uigetdir; %Change to where you stored the .mat files
%newfilename='Jason_data'; %Change to the desired name you want save with.
[newfilename,saveloc]=uiputfile('*.mat'); %Change to the location where you want to save to
listing=dir([fileloc '\*.mat']);

if exist([newfilename '.mat'],'file')
    load([newfilename '.mat']);
    g=size(s,1);
else
    s=struct('Name',{},'FireRate',[],'Pulses',{},'times',{},'clusters',{},'waveforms',{},'Stim',{},'Intensity',{},'Area',{},'Depth',{});
    g=0;
end
for k=g+1:g+size(listing,1);
    file=listing(k).name(1:end-4);
    load([fileloc '\' file '.mat'],'-regexp', '_Ch');
    %I'm going to assume that all waveforms are in Ch4, Ch2 or Ch3
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
    if exist(['V' file '_Ch32']) & eval(['size(V' file '_Ch32.times,1)>0'])
        eval(['time_trig=V' file '_Ch32.times']);
        time_trig=time_trig+0.00282*ones(length(time_trig),1);
        s(k).Pulses=1000*time_trig;
    elseif sum(data.codes(:,1)==9)>0
        s(k).Pulses=1000*data.times(data.codes(:,1)==9);
    else
        figure(1)
        for b=max(data.codes(:,1)):-1:1
            subplot(2,5,double(b))
            plot(data.values(data.codes(:,1)==b,:)')
            title(sprintf('Code is %d',b))
        end
        for m=max(data.codes(:,1)):-1:1
            if sum(data.codes(:,1)==m)>0
                figure(2)
                plot(data.values(data.codes(:,1)==m,:)')
                check=input('Is this the TMS pulse?\n');
                if check==1
                    s(k).Pulses=1000*data.times(data.codes(:,1)==m);
                    sprintf('Exit \n')
                    clf(1)
                    clf(2)
                    break
                else
                    clf(2)
                end
                sprintf('Still Here \n')
            end
        end
        clf
    end
    clear -regexp '_Ch'
end

save([saveloc newfilename],'s')