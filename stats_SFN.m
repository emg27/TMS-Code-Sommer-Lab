%%Create the statistics for the SFN poster and the TMS paper. The current
%stats we plan to run: 
%
% 1) To test if there is a significant difference subthres for stim vs sham
% 2) To test if there is a significant difference suprathres for stim/sham
% 3) Compare the differences between the stim responses (sub vs supra)
% 4) Stats for the period of inhibition (bootstrapping)
% 5) Stats for the period of excitation (bootstrapping)

%%Significant difference between subthreshold for stim and sham; take the
%%difference between the two wave forms and then take a paired t-test

%Stim=avgSham(1:9,(tbase+1)-50:(tbase+1)+ta);
pvalue=.05/1000;

crop = (tbase+1)-50:(tbase+1)+ta;
compcrop = (tbase+1):(tbase+1)+ta;
Trunk=3+gauss_size:size(shamps,2)-gauss_size;

wStimSub=stimps(stimps(:,2)<=50,Trunk);
wStimSup=stimps(stimps(:,2)>50,Trunk);
wShamSub=shamps(shamps(:,2)<=50,Trunk);
wShamSup=shamps(shamps(:,2)>50,Trunk);

%Create the difference plots
diffSub=mean(wStimSub)-mean(wShamSub);
diffSup=mean(wStimSup)-mean(wShamSup);
diffStim=mean(wStimSup)-mean(wStimSub);
diffSham=mean(wShamSup)-mean(wShamSub);


%1)
[Hsub,Psub]=ttest2(wStimSub(:,compcrop),wShamSub(:,compcrop));

%2) 
[Hsup,Psup]=ttest2(wStimSup(:,compcrop),wShamSup(:,compcrop));

%3)
[HcompSt,PcompSt]=ttest2(wStimSup(:,compcrop),wStimSub(:,compcrop));

%4)
[HcompSh,PcompSh]=ttest2(wShamSup(:,compcrop),wShamSub(:,compcrop));

%Plot the difference
figure
subplot(2,2,1)
plot(-50:ta,diffSub(crop),'g', -50:ta, mean(wStimSub(:,crop)), 'b',...
    -50:ta, mean(wShamSub(:,crop)) ,'r',...
    compcrop(Psub<pvalue)-(tbase+1),-0.4*ones(size(compcrop(Psub<pvalue))),'k*')
legend('Difference','Stim','Sham')
title('Subthreshold Difference')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

subplot(2,2,3)
plot(-50:ta,diffSup(crop),'g', -50:ta, mean(wStimSup(:,crop)), 'b',...
    -50:ta, mean(wShamSup(:,crop)) ,'r',...
    compcrop(Psup<pvalue)-(tbase+1),-0.4*ones(size(compcrop(Psup<pvalue))),'k*')
legend('Difference','Stim','Sham')
title('Suprathreshold Difference')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

subplot(2,2,2)
plot(-50:ta,diffStim(crop),'g', -50:ta, mean(wStimSup(:,crop)), 'b',...
    -50:ta, mean(wStimSub(:,crop)) ,'c',...
    compcrop(PcompSt<pvalue)-(tbase+1),-0.4*ones(size(compcrop(PcompSt<pvalue))),'k*')
legend('Difference','Suprathreshold','Subthreshold')
title('Sub vs Suprathreshold Difference for Stim')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

subplot(2,2,4)
plot(-50:ta,diffSham(crop),'g', -50:ta, mean(wShamSup(:,crop)), 'b',...
    -50:ta, mean(wShamSub(:,crop)) ,'c',...
    compcrop(PcompSh<pvalue)-(tbase+1),-0.4*ones(size(compcrop(PcompSh<pvalue))),'k*')
legend('Difference','Suprathreshold','Subthreshold')
title('Sub vs Suprathreshold Difference for Sham')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

%%Create a sliding window average of 50ms of data
n=0;
window=20; %ms
startpt=tbase+1;
while (startpt+window)<(tbase+ta+1)
n=n+1;
WindTime(n)=n*window;
StSubWin(:,n)=mean(wStimSub(:,startpt:startpt+window)')';
StSupWin(:,n)=mean(wStimSup(:,startpt:startpt+window)')';
ShSubWin(:,n)=mean(wShamSub(:,startpt:startpt+window)')';
ShSupWin(:,n)=mean(wShamSup(:,startpt:startpt+window)')';
startpt=tbase+1+ceil(n*window/2);

end

%Run TTest over the window averaged data
%1)
[HsubW,PsubW]=ttest2(StSubWin,ShSubWin);
%2) 
[HsupW,PsupW]=ttest2(StSupWin,ShSupWin);
%3)
[HcompStW,PcompStW]=ttest2(StSupWin,StSubWin);
%4)
[HcompShW,PcompShW]=ttest2(ShSupWin,ShSubWin);

%Plot the difference
figure
subplot(2,2,1)
plot(-50:ta,diffSub(crop),'g', -50:ta, mean(wStimSub(:,crop)), 'b',...
    -50:ta, mean(wShamSub(:,crop)) ,'r',...
    WindTime(PsubW<pvalue),-0.4*ones(size(WindTime(PsubW<pvalue))),'k*')
legend('Difference','Stim','Sham')
title(['Subthreshold Difference with Sliding Window=' num2str(window) 'ms'])
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

subplot(2,2,3)
plot(-50:ta,diffSup(crop),'g', -50:ta, mean(wStimSup(:,crop)), 'b',...
    -50:ta, mean(wShamSup(:,crop)) ,'r',...
    WindTime(PsupW<pvalue),-0.4*ones(size(WindTime(PsupW<pvalue))),'k*')
legend('Difference','Stim','Sham')
title('Suprathreshold Difference')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

subplot(2,2,2)
plot(-50:ta,diffStim(crop),'g', -50:ta, mean(wStimSup(:,crop)), 'b',...
    -50:ta, mean(wStimSub(:,crop)) ,'c',...
    WindTime(PcompStW<pvalue),-0.4*ones(size(WindTime(PcompStW<pvalue))),'k*')
legend('Difference','Suprathreshold','Subthreshold')
title('Sub vs Suprathreshold Difference for Stim')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

subplot(2,2,4)
plot(-50:ta,diffSham(crop),'g', -50:ta, mean(wShamSup(:,crop)), 'b',...
    -50:ta, mean(wShamSub(:,crop)) ,'c',...
    WindTime(PcompShW<pvalue),-0.4*ones(size(WindTime(PcompShW<pvalue))),'k*')
legend('Difference','Suprathreshold','Subthreshold')
title('Sub vs Suprathreshold Difference for Sham')
xlabel('Time (ms)')
axis([-50 ta -0.4 0.4])

%%Extra assignment: plot the mean firing rate of dur(T) vs intensity and
%%dur(T) length.


% Width=5:10:300; %increase averaging window by 10 
% 
% [inten,width]=meashgrid(0:10:90,Width);
% Fire=0.*inten(:);
% 
% for n=1:length(Fire)
% end
