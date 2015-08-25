%% Find the peak value immediately following the TMS pulse.
close all
thres=0.01;
thresnorm=0.4;
bin_size=6;
endpt=101; 

% Calculate the values for all the peaks following the TMS pulse
% timeVal: Times at which peaks occur
% peakVal: values of peaks when they occur
peakVal=[];
peakVal=[];
allpostTMS=allptsh(:,3+gauss_size+tbase+1:end-gauss_size);
normpostTMS=normptsh(:,3+gauss_size+tbase+1:end-gauss_size);
check=nan(n,12);
ct=check;
for n=1:size(allptsh,1)
    [peak,time]=findpeaks(allpostTMS(n,1:101))%,'MinPeakHeight',thres);
    %[peak,time]=findpeaks(allptsh(n,3+tbase+gauss_size:end-gauss_size),'MinPeakHeight',thres);
    [pnorm,tnorm]=findpeaks(normpostTMS(n,1:101),'MinPeakHeight',thresnorm);
    %[pnorm,tnorm]=findpeaks(normptsh(n,3+tbase+gauss_size:end-gauss_size),'MinPeakHeight',thresnorm);
    check(n,1:length(pnorm))=pnorm;
    ct(n,1:length(pnorm))=tnorm;
    if size(time,2)>0
        peakVal(n,1)=peak(1);%max(peak);%peak(1)
        %xPeak=find(allpostTMS(n,1:101)==max(peak));
        timeVal(n,1)=time(1);%xPeak(1);%time(1);
    else
        peakVal(n,1)=nan;
        timeVal(n,1)=nan;
    end
    
    if size(tnorm,2)>0
        peakVal(n,2)=pnorm(1);%max(pnorm);%pnorm(1);
        xPnorm=find(normpostTMS(n,1:101)==max(pnorm));
        timeVal(n,2)=tnorm(1);%xPnorm(1);%tnorm(1);
    else
        peakVal(n,2)=nan;
        timeVal(n,2)=nan;
    end
end

%Determine the histogram distribution for the stim files
[StBin,timeSt]=hist(timeVal(pSt,1),0:bin_size:endpt);
[StnormBin,tnormSt]=hist(timeVal(pSt,2),0:bin_size:endpt);
%Determine the histogram distribution for the sham files
[ShBin,timeSh]=hist(timeVal(pSh,1),0:bin_size:endpt);
[ShnormBin,tnormSh]=hist(timeVal(pSh,2),0:bin_size:endpt);

figure
subplot(2,1,1)
% plot(timeSt,StBin,'bo-',timeSh,ShBin,'go-')
% legend('Stim','Sham')
% subplot(2,1,2)
plot(tnormSt,StnormBin,'bo-',tnormSh,ShnormBin,'go-')
legend('Stim Normalized','Sham Normalized')
xlim([0 200])
subplot(2,1,2)
plot(tnormSt,StnormBin/sum(StnormBin),'bo-',...
    tnormSh,ShnormBin/sum(ShnormBin),'go-')
legend('Stim Normalized','Sham Normalized')
xlim([0 200])

% intenSt=figure;
intenN=figure;
% intenSt2=figure;
intenN2=figure;

for n=1:9
    posSt=find(normptsh(:,2)<=n*10 & normptsh(:,2)>10*(n-1) & normptsh(:,1)==1);
    posSh=find(normptsh(:,2)<=n*10 & normptsh(:,2)>10*(n-1) & normptsh(:,1)==0);
    peakmSt(n,1:2)=nanmean(peakVal(posSt,:));
    peakmSh(n,1:2)=nanmean(peakVal(posSh,:));
    timeMSt(n,1:2)=nanmean(timeVal(posSt,:));
    timeMSh(n,1:2)=nanmean(timeVal(posSh,:));
    %Determine the histogram distribution for the stim files
    [StBin(n,:),timeSt(n,:)]=hist(timeVal(posSt,1),0:bin_size:endpt);
    [StnormBin(n,:),tnormSt(n,:)]=hist(timeVal(posSt,2),0:bin_size:endpt);
    %Determine the histogram distribution for the sham files
    [ShBin(n,:),timeSh(n,:)]=hist(timeVal(posSh,1),0:bin_size:endpt);
    [ShnormBin(n,:),tnormSh(n,:)]=hist(timeVal(posSh,2),0:bin_size:endpt);
%     figure(intenSt)
%     subplot(3,3,n)
%     plot(timeSt(n,:),StBin(n,:),'bo-',timeSh(n,:),ShBin(n,:),'go-')
%     title([num2str(n*10) '% Stim'])
%     xlim([0 150])
    figure(intenN)
    subplot(3,3,n)
    plot(tnormSt(n,:),StnormBin(n,:),'bo-',tnormSh(n,:),ShnormBin(n,:),'go-')
    title([num2str(n*10) '% Stim Normalized'])
    xlim([0 100])
%     figure(intenSt2)
%     subplot(3,3,n)
%     plot(timeSt(n,:),StBin(n,:)/sum(StBin(n,:)),'bo-',...
%         timeSh(n,:),ShBin(n,:)/sum(ShBin(n,:)),'go-')
%     title([num2str(n*10) '% Stim Percentage'])
    xlim([0 150])
    figure(intenN2)
    subplot(3,3,n)
    plot(tnormSt(n,:),StnormBin(n,:)/sum(StnormBin(n,:)),...
        'bo-',tnormSh(n,:),ShnormBin(n,:)/sum(ShnormBin(n,:)),'go-')
    title([num2str(n*10) '% Stim Normalized Percentage'])
    xlim([0 150])
end
figure;
plot(10:10:90,timeMSt(:,2),'ro',10:10:90,timeMSh(:,2),'bo')
% %% ANOVAN for binned values across: time, Intensity, and type
% % Anova Parameters
% gInt = 10:10:90;
% gTimes = 0:6:200;
% gType = 0:1;
% 
% % 3 dim matrix with all values.
% %   1: Intensity
% %   2: Times
% %   3: Stim/Sham
% dataAnova(:,:,1) = StnormBin(:,1:length(gTimes));
% dataAnova(:,:,2) = ShnormBin(:,1:length(gTimes));
% 
% % Vectorize Anova
% dataAnovaVector = [];
% dataAnovaVectorInd = 1;
% gIntVector = [];
% gTimesVector = [];
% gTypeVector = [];
% 
% % Vectorize
% % i: Stimulation type
% for i = 1:size(dataAnova,3)
%     % j: Time information
%     for j = 1:size(dataAnova,2)
%         % k: Intensity
%         for k = 1:size(dataAnova,1)
%             % Vector:
%             % stim sham stim sham stim sham stim sham stim sham stim sham
%             % tim1 tim1 tim2 tim2 tim3 tim3 tim4 tim4 tim5 tim5 tim6 tim6
%             % int1 int1 int1 int1 int1 int1 int1 int1 int1 int1 int1 int1
%             dataAnovaVector(dataAnovaVectorInd) = dataAnova(k,j,i);
%             gIntVector(dataAnovaVectorInd) = gInt(k);
%             gTimesVector(dataAnovaVectorInd) = gTimes(j);
%             gTypeVector(dataAnovaVectorInd) = gType(i);
%             dataAnovaVectorInd = dataAnovaVectorInd + 1;
%         end
%     end
% end
% 
% % Perform Anova
% [A, B, stats] = anovan(dataAnovaVector,...
%     {gTypeVector', gTimesVector', gIntVector'})
% 
% %% Anova for peak times across: Intensity, stim/sham
% yTimes = timeVal(:,1);
% gType = allptsh(:,1);
% gInt = allptsh(:,2);
% 
% % Remove NAN for timeVals
% temp=yTimes;
% yTimes=yTimes(~isnan(temp)); gType=gType(~isnan(temp)); gInt=gInt(~isnan(temp));
% % Remove times > 200 ms for timeVals
% temp=yTimes;
% yTimes=yTimes(temp < 100); gType=gType(temp < 100); gInt=gInt(temp < 100);
% % Remove NAN for Type
% temp=gType;
% yTimes=yTimes(~isnan(temp)); gType=gType(~isnan(temp)); gInt=gInt(~isnan(temp));
% % Remove unknown stimulus conditions
% temp=gType;
% yTimes=yTimes(temp ~= 3); gType=gType(temp ~= 3); gInt=gInt(temp ~= 3);
% % Remove NAN for Intensity conditionsclear
% 
% temp=gInt;
% yTimes=yTimes(~isnan(temp)); gType=gType(~isnan(temp)); gInt=gInt(~isnan(temp));
% 
% % Perform ANOVA
% [A, B, stats] = anovan(yTimes, {gType, gInt}, 'model', 'interaction');
% 
% %% Anova for Stim for peak times across: Intensity
% % Remove sham conditions
% temp=gType;
% yTimes=yTimes(temp ~= 0); gType=gType(temp ~= 0); gInt=gInt(temp ~= 0);
% [A, B, stats] = anovan(yTimes, {gInt});