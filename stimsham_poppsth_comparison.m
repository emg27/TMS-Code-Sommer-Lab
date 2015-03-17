close all
tb=1000;
ta=1000;
gauss=5;
boy=ones(375,1)*median(widthall2);
temp=widthall2<boy;
for k=1:9
%intensity=50

wdthSt=[];
wdthSh=[];
mfr_Stim=[];
mfr_Sham=[];
inten_stim=[];
inten_sham=[];
check=[];
cellSt=0;
cellSh=0;
aveWidth=mean(temp,2);
includeCells=find(stiminten2(:,2)<=k*10 & stiminten2(:,2)>(k-1)*10);
avWdth=aveWidth(includeCells);
for n=1:size(includeCells,1)
    date=stiminten2(includeCells(n),3);
    cluster=stiminten2(includeCells(n),4);
    stim=stiminten2(includeCells(n),1);
    
    if stim==1
        color='r';
    elseif stim==0
        color='b';
    else
        continue
    end
    figure(1)
    hold on
    clust=find(s(date).clusters==cluster);
    [spk_d,trl_fr,bin_t,baseline,m_fr,~]=psth1block(s(date).Pulses,tb,ta,1000*s(date).times(clust),gauss,0);
    set(spk_d,'Color',color)
    if stim==1
        if size(m_fr,2)>1
            mfr_Stim=[mfr_Stim; 1000*m_fr];
            wdthSt=[wdthSt; avWdth(n)];
            cellSt=cellSt+1;
        end
    elseif stim==0
        if size(m_fr,2)>1
            mfr_Sham=[mfr_Sham; 1000*m_fr];
            cellSh=cellSh+1;
            wdthSh=[wdthSh; avWdth(n)];
        end
    end
    check=[check; s(date).Stim{1}];
end

%%All together
figure
subplot(2,1,1)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(mfr_Stim),'Color','r')
title(num2str(cellSt))
subplot(2,1,2)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(mfr_Sham),'Color','b')
title(num2str(cellSh))
xlabel([num2str(n*10) '%'])

%%Broken down
inhibSt=mfr_Stim(wdthSt<.4,:);
inhibSh=mfr_Sham(wdthSh<.4,:);
excitSt=mfr_Stim(wdthSt>=.7,:);
excitSh=mfr_Sham(wdthSh>=.7,:);
orSt=mfr_Stim(wdthSt>=.4 & wdthSt<.7,:);
orSh=mfr_Sham(wdthSh>=.4 & wdthSh<.7,:);

%mean
inhibStm=mean(inhibSt);
inhibShm=mean(inhibSh);
excitStm=mean(excitSt);
excitShm=mean(excitSh);
orStm=mean(orSt);
orShm=mean(orSh);

%Standard Error
inhibStse=std(inhibSt)./sqrt(size(inhibSt,1));
inhibShse=std(inhibSh)./sqrt(size(inhibSh,1));
excitStse=std(excitSt)./sqrt(size(excitSt,1));
excitShse=std(excitSh)./sqrt(size(excitSh,1));
orStse=std(orSt)./sqrt(size(orSt,1));
orShse=std(orSh)./sqrt(size(orSh,1));


figure
subplot(3,2,1)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(inhibSt),'r')
title(num2str(length(wdthSt(wdthSt<.4,:))))

subplot(3,2,2)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(mfr_Sham(wdthSh<.4,:)),'b')
title(num2str(length(wdthSh(wdthSh<.4,:))))

subplot(3,2,3)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(orSt),'r')
title(num2str(length(wdthSt(wdthSt>=.4 & wdthSt<.7))))

subplot(3,2,4)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(orSh),'b')
title(num2str(length(wdthSh(wdthSh>=.4 & wdthSh<.7))))

subplot(3,2,5)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(excitSt),'r')
title(num2str(length(wdthSt(wdthSt>=.7,:))))

subplot(3,2,6)
plot(bin_t,0.*bin_t,'k--',bin_t,mean(excitSh),'b')
title(num2str(length(wdthSh(wdthSh>=.7,:))))
xlabel([num2str(k*10) '%'])

% figure
% subplot(3,2,1)
% 
% hold on
% plot_variance(bin_t,inhibStm-inhibStse,inhibStm+inhibStse,'r');
% plot(bin_t,0.*bin_t,'k--',bin_t,mean(inhibSt),'k')
% title(num2str(length(wdthSt(wdthSt<.4,:))))
% subplot(3,2,2)
% 
% hold on
% plot_variance(bin_t,inhibShm-inhibShse,inhibShm+inhibShse,'b');
% plot(bin_t,0.*bin_t,'k--',bin_t,mean(mfr_Sham(wdthSh<.4,:)),'k')
% title(num2str(length(wdthSh(wdthSh<.4,:))))
% 
% subplot(3,2,3)
% 
% hold on
% plot_variance(bin_t,excitStm-excitStse,excitStm+excitStse,'r');
% plot(bin_t,0.*bin_t,'k--',bin_t,mean(orSt),'k')
% title(num2str(length(wdthSt(wdthSt>=.4 & wdthSt<.7))))
% subplot(3,2,4)
% 
% hold on
% plot_variance(bin_t,excitShm-excitShse,excitShm+excitShse,'b');
% plot(bin_t,0.*bin_t,'k--',bin_t,mean(orSh),'k')
% title(num2str(length(wdthSh(wdthSh>=.4 & wdthSh<.7))))
% 
% subplot(3,2,5)
% 
% hold on
% plot_variance(bin_t,orStm-orStse,orStm+orStse,'r');
% plot(bin_t,0.*bin_t,'k--',bin_t,mean(excitSt),'k')
% title(num2str(length(wdthSt(wdthSt>=.7,:))))
% subplot(3,2,6)
% 
% hold on
% plot_variance(bin_t,orShm-orShse,orShm+orShse,'b');
% plot(bin_t,0.*bin_t,'k--',bin_t,mean(excitSh),'k')
% title(num2str(length(wdthSh(wdthSh>=.7,:))))
% xlabel([num2str(k*10) '%'])
end