%%Create the Poisson Spike Train Analysis to determine how improbable it is
%%that the number of action potentials within a specific time interval is a
%%chance occurance. We will be comparing actual spikes to the number of
%%predicted spikes from a Poisson distributed for a given time interval.

%Clear and set up the workspace
clear
[filename, pathname]=uigetfile('*.mat'); %Load the structure
load([pathname filename])
uiload %load the TMS neuron counts
close all

%Determine the time range of each trial:
%I will initially set each trial period as 0.5 seconds before and after each
%TMS pulse.
ta=300; %Time after the TMS pulse in ms
tb=100; %Time before the TMS pulse in ms
allburst=[];
Name=cell(0);
spkburstavg=nan(size(s,2),9);
blackouts=nan(size(s,2),1);
stimsham=nan(size(s,2),1);
%Set equations
rho=@(lambda,spknum) 1-cdf('Poisson',spknum,lambda);
SI=@(lambda,spknum) -log(rho(lambda,spknum));
count=0;
for k=1:size(s,2)
    block=s(k);
    cluster=unique(block.clusters);
    
    for neur=1:length(cluster)
        if cluster(neur)==0 | cluster(neur)==9
            continue
        end
    neuron=cluster(neur);
    Pulses=block.Pulses;
    close all
    figure
    if length(block.Pulses)>0 & size(block.Intensity{1})>0
        clusterpos=find(block.clusters==neuron);
        if size(clusterpos)<1 | ~strcmp(s(k).BrainArea{1}(1),'M')
            continue
        end
    else
        continue
    end
    burst=nan(length(Pulses),9);
    count=count+1;
    for n=1:length(Pulses)
        
        %Determine if the block is Sham, Stim, or undefined
        if strfind(s(k).Stim{:},'Sh')==1 %Used Sh in case of typo for Sham
            StSh=0;
        elseif strfind(s(k).Stim{:},'St')==1 %Used St in case of typo for Stim
            StSh=1;
        else
            StSh=-1;
        end
        %             burst(n,4:6)=[str2num(inten{1})/10 strcmp(s(k).Stim{1},'Stim') k];
        Name=[Name; s(k).Name];
        inten=block.Intensity(1);
        NeurNum=AllInfo(AllInfo(:,2)==k,7);
        burst(n,4:7)=[str2num(inten{1}) StSh k NeurNum];
        
        if n>length(Pulses)
            continue
        end
        %Find spike times
        [position timepts]=Raster(Pulses(n),tb,ta,1000*block.times(clusterpos),n);
        
        %close
        
        %Set index spike as the first spike after the TMS pulse
        spkinTrial=block.times(clusterpos(timepts))-Pulses(n)/1000;%s
        spk_aft_TMS=spkinTrial(spkinTrial>=0);
        index=1;
        
        %Calculate r
        r=1000*length(timepts)/(ta+tb); %seconds
%         r=1000*length(spkinTrial(spkinTrial<0))/ta;

        if length(spk_aft_TMS)>0
            burst(n,8)=spk_aft_TMS(1);
            if size(str2num(s(k).Blackout{1}),2)>0
                burst(n,9)=str2num(s(k).Blackout{1});
            else
                burst(n,9)=nan;
            end
        end
        if length(spk_aft_TMS)<=2
            continue
        end
        %Find the two consecutive spikes with rate>=r
        isiIndex=diff(spk_aft_TMS(index:end));
        isi_GrEq_r_pos=find((1./isiIndex)>r);
        if size(isi_GrEq_r_pos,1)>0
            %Tbegin=0;
            TBrstEnd=isiIndex(isi_GrEq_r_pos(1));
            TBrstBeg=diff(spk_aft_TMS([index end]));
            spkNum=1;
            [Tall,rhoall]=deal(nan(length(isiIndex),1));
            SIall=nan(length(isiIndex),2);
            
            for isi=1:length(isiIndex)
                if isi>=(isi_GrEq_r_pos(1)+1)
                    %Find the end of the burst period
                    TBrstEnd=TBrstEnd+isiIndex(isi);
                    Tall(isi)=TBrstEnd;
                    rhoall(isi)=rho(r*TBrstEnd,spkNum);
                    SIall(isi,1)=SI(r*TBrstEnd,spkNum);
                    spkNum=spkNum+1;
                end
            end
            [valMax posBrstEnd]=max(SIall(:,1));
            posBrstEnd=posBrstEnd+1;
            
            for neuron=1:posBrstEnd;
            %Find the begining of the burst period
            SIall(neuron,2)=SI(r*TBrstEnd,length(spk_aft_TMS(neuron:end)));
            end
            
            [valMax posMax]=max(SIall);
            burst(n,1:3)=[spk_aft_TMS(isi_GrEq_r_pos(1))...
                spk_aft_TMS([posMax(1) posMax(2)])'];
            
           
        end
        
        hold on
%         plot(burst(n,3)*1000,n,'k*',burst(n,2)*1000,n,'bo')
        plot(burst(n,3)*1000,n,'ro','MarkerFaceColor','r')
        
    end
    if size(str2num(s(k).Blackout{1}),2)>0
        rectangle('Position',[0 0 str2num(s(k).Blackout{1}) n],...
            'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75])
        blackouts(k)=str2num(s(k).Blackout{1});
    end
    axis([-tb ta 0 n])
    xlabel('Time (ms)'),ylabel('Pulse Number')
        title([s(k).Name ' ' s(k).Stim{1} ' ' s(k).Intensity{1} ...
            '% Neuron ' num2str(NeurNum) ' out of ' num2str(max(AllInfo(:,6)))])
        %print(s(k).Name(1:8) '_' s(k).Stim{1}...
         %   num2str(burst(n,4)*10) '_Neuron' num2str(neuron)])
%          print(gcf, '-dpng',['Neuron' num2str(count)])
        %close 
    if strcmp(s(k).Stim{1},'Stim')
        stimsham(k)=1;
    elseif strcmp(s(k).Stim{1},'Sham')
        stimsham(k)=0;
    end
    if size(burst)>0
        allburst=[allburst; burst];
        spkburstavg(k,:)=nanmean(burst(:,1:9));
    end
    end
end

figure
subplot(3,1,1)
hist(spkburstavg(~isnan(spkburstavg(:,1)),3),25)
title('Calculated Burst Onset time for all Stim/Sham Intensities')
xlabel('Time (sec)'),ylabel('Block Count')
ylim([0 300])

subplot(3,1,2)
hist(spkburstavg(~isnan(spkburstavg(:,1)),2),25)
title('Calculated Burst Offset time for all Stim/Sham Intensities')
xlabel('Time (sec)'),ylabel('Block Count')
ylim([0 300])

subplot(3,1,3)
hist(spkburstavg(~isnan(spkburstavg(:,1)),1),25)
title('Second Spike in the Burst"s Onset time for all Stim/Sham Intensities')
xlabel('Time (sec)'),ylabel('Block Count')
ylim([0 300])

binnum=0:5:300;
stimOn=figure;shamOn=figure;
stimOnVar=figure;shamOnVar=figure
for n=1:9
    figure(stimOn)
    subplot(3,3,n)
    posStim=find(allburst(:,4)==n*10 & allburst(:,5)==1);
    hist(1000*allburst(posStim,3),binnum)
    title(['Calculated Burst Onset time for Stim ' num2str(n*10) '%'])
    xlabel('Time (ms)'),ylabel('Block Count')
    axis([0 300 0 350])
%     figure(stimOnVar)
%     subplot(3,3,n)
%     hist(nanstd(spkburstavg(posStim,3)),binnum)
%     title(['Standard Deviation of Calculated Burst Onset time for Stim ' num2str(n*10) '%'])
%     xlabel('Time (sec)'),ylabel('Block Count')

    figure(shamOn)
    subplot(3,3,n)
    posSham=find(allburst(:,4)==n*10 & allburst(:,5)==0);
    hist(1000*allburst(posSham,3),binnum)
    title(['Calculated Burst Onset time for Sham ' num2str(n*10) '%'])
    xlabel('Time (ms)'),ylabel('Block Count')
    axis([0 300 0 150])
    
    figure(stimOnVar)
    hold on
    plot(n*10,nanstd(1000*allburst(posStim,3)),'ko')
    plot(n*10,nanstd(1000*allburst(posSham,3)),'ro')
    title('Standard Deviation of Calculated Burst Onset for Stim and Sham')
    xlabel('Intensity (\%)');ylabel('Time (ms)')
    legend('Stim','Sham')
    
%     figure(shamOnVar)
%     subplot(3,3,n)
%     hist(nanstd(spkburstavg(posSham,3)),binnum)
%     title(['Standard Deviation of Calculated Burst Onset time for Sham ' num2str(n*10) '%'])
%     xlabel('Time (sec)'),ylabel('Block Count')

end

PositionInStructure=allburst(:,6); NeuronNumber=allburst(:,7); 
FirstSpike=1000*allburst(:,8); SecondSpike=1000*allburst(:,1); 
Onset=1000*allburst(:,3); Offset=1000*allburst(:,2); 
Intensity=allburst(:,4); StimSham=allburst(:,5);
Blackout=allburst(:,9);
T=table(PositionInStructure,Name,NeuronNumber,Intensity,StimSham,Onset,Offset,FirstSpike,SecondSpike, Blackout);
filename='PoissonSpikeTrainDataNew_wBlackout.xlsx';
writetable(T,filename,'Sheet',1,'Range','B1')