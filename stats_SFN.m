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

crop = (tbase+1)-50:(tbase+1)+ta;
compcrop = (tbase+1):(tbase+1)+ta;
wStimSub=avgStim_sub;
wStimSup=avgStim_supra;
wShamSub=avgSham_sub;
wShamSup=avgSham_supra;

%1)
diffSub=wStimSub-wShamSub
%[Hsub(1),Psub(1)]=ttest(diffSub(compcrop),zeros(size(diffSub(compcrop))));
[Hsub,Psub]=ttest2(wStimSub(compcrop),wShamSub(compcrop));
%2) 
diffSup=wStimSup-wShamSup
%[Hsup(2),Psup(2)]=ttest(diffSup(compcrop),zeros(size(diffSub(compcrop))));
[Hsup,Psup]=ttest2(wStimSup(compcrop),wShamSup(compcrop));

%3)
diffStim=wStimSup-wStimSub
%[Hcomp(3),Pcomp(3)]=ttest(diffStim(compcrop),zeros(size(diffStim(compcrop)))); 
[HcompSt,PcompSt]=ttest2(wStimSup(compcrop),wStimSub(compcrop))

%4)
diffSham=wShamSup-wShamSub
%[Hcomp(4),Pcomp(4)]=ttest(diffSham(compcrop),zeros(size(diffSham(compcrop)))); 
[HcompSh,PcompSh]=ttest2(wShamSup(compcrop),wShamSub(compcrop))

%Plot the difference

figure
subplot(2,2,1)
plot(-50:ta,diffSub(crop),'g', -50:ta, wStimSub(crop), 'b', -50:ta, wShamSub(crop) ,'r')
legend('Difference','Stim','Sham')
title('Subthreshold Difference')
xlabel(['p=' num2str(Psub)])
xlim([-50 ta])

subplot(2,2,3)
plot(-50:ta,diffSup(crop),'g', -50:ta, wStimSup(crop), 'b', -50:ta, wShamSup(crop) ,'r')
legend('Difference','Stim','Sham')
title('Suprathreshold Difference')
xlabel(['p=' num2str(Psup)])
xlim([-50 ta])

subplot(2,2,2)
plot(-50:ta,diffStim(crop),'g', -50:ta, wStimSup(crop), 'b', -50:ta, wStimSub(crop) ,'c')
legend('Difference','Suprathreshold','Subthreshold')
title('Sub vs Suprathreshold Difference for Stim')
xlabel(['p=' num2str(PcompSt)])
xlim([-50 ta])

subplot(2,2,4)
plot(-50:ta,diffSham(crop),'g', -50:ta, wShamSup(crop), 'b', -50:ta, wShamSub(crop) ,'c')
legend('Difference','Suprathreshold','Subthreshold')
title('Sub vs Suprathreshold Difference for Sham')
xlabel(['p=' num2str(PcompSh)])
xlim([-50 ta])
%%Extra assignment: plot the mean firing rate of dur(T) vs intensity and
%%dur(T) length.


% Width=5:10:300; %increase averaging window by 10 
% 
% [inten,width]=meashgrid(0:10:90,Width);
% Fire=0.*inten(:);
% 
% for n=1:length(Fire)
% end
