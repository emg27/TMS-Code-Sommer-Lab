clear
[filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
load([pathname filename])
close all

widthall=[];
cnta=0
baseline=[];
spks=[];
guess=[];
Erinn=[];
spksvar=[];
i=1;
a=[];
stiminten=double([]);
gue=figure
for n=1:size(s,2)
    if length(s(n).Pulses)>0 
        for g=1:max(s(n).clusters)
            if size(s(n).width,1)==0
            % if s(n).Good(g)~=1 | size(s(n).width,1)==0
                continue
            end
                try2mak=[];
                tempEri=[n double(g)];
                Erinn=[Erinn; tempEri];
                clust=find(s(n).clusters==g);
                if length(clust)>70 && length(clust)~=length(s(n).Pulses)
                plot(s(n).waveforms(clust,:)')
                if length(clust)>1
                    sptemp=mean(s(n).waveforms(clust,:));
                else
                    sptemp=s(n).waveforms(clust,:);
                end
                for try1=1:length(clust)
                        try2mak(try1,:)=s(n).waveforms(clust(try1),:)/max(s(n).waveforms(clust(try1),:));
                end
                spvar=var(try2mak);
                time=(0:1:length(sptemp)-1)/50;
                spikes=nan(1,300);
                spikesV=nan(1,300);
                sptemp=sptemp/max(abs(sptemp));
                a=[a; mean(try2mak)' sptemp'];
                spvar=spvar;
                pos=find(sptemp==min(sptemp));
                spikes(150-pos:150-pos-1+length(sptemp))=sptemp;
                spikesV(150-pos:150-pos-1+length(sptemp))=spvar;
%                 if size(widthall,2)<=0
%                     widthall=s(n).width{g};
%                     baseline=s(n).Base(g);
%                     spks=spikes;
%                    % guess=thought;
%                 end
                if size(s(n).width{g},1)~=1
                    s(n).width{g}=s(n).width{g}';
                end
            spks=[spks; spikes]; %[spks; size(s(n).waveforms(clust,:))];
            spksvar=[spksvar; spikesV];
            %guess=[guess; thought];
            i=size(spks,1);
            %pause
            widthall=[widthall; s(n).width{g}];
            baseline=[baseline; 1000*s(n).Base(g)];
            temp=nan(length(baseline),9);
%             for k=1:9
%                 temp(i,k)=(slope(k,1)*widthall(i,k)+slope(k,2))>baseline(i);
%             end
            guess2=mean(temp(i,:));
            %l(i)=guess2;
%             if length(find(axon==size(spks,1)))>0
%                 sprintf('%s C%d Spks=%d',s(n).Name,g,size(spks,1)) 
%             end
%            if find(axon==i)>0
%                 s(n).type(g)='A';
%                 cnta=cnta+1;
            if guess2<=.3 & guess2>=0
                s(n).type(g)='I';
            elseif guess2>.7
                s(n).type(g)='E';
            elseif guess2>.3 && guess2<=.7
                s(n).type(g)='U'; %u is for unknown
            end
            if strcmp(s(n).Stim{1},'Sham')==1
                stim=0;
            elseif strcmp(s(n).Stim{1},'Stim')==1
                stim=1;
            else
                stim=3;
            end
            
            stiminten=[stiminten; stim s(n).Intensity(1) n double(g)];
%            else 
%                s(n).type(g)='_';
%            end
                end
        end
    end
end

baseline=baseline(~any(isnan(widthall), 2), :);
widthall= widthall(~any(isnan(widthall), 2), :);
stiminten=stiminten(~any(isnan(widthall), 2), :);

remove8=find(widthall(:,8)~=0 & baseline>0);
widthall2=widthall(remove8,:);
stiminten2=stiminten(remove8,:);
baseline2=baseline(remove8,:);
temp=nan(length(baseline2),10);
for k=1:9
%     [width, time]=hist(widthall2(:,k),30);
    [idx,C]=kmeans(widthall2(:,k),2,'Distance','cityblock');
    temp(:,k)=idx;
end
[idx,C]=kmeans(log(baseline2),2,'Distance','cityblock');
temp(:,10)=idx;
temp(:,10)=abs(temp(:,10)-3);
guess2=mean(temp,2)-1;
%guess2(axon)=-1*(guess2(axon)+ones(length(axon),1));

for k=1:9
    figure(1)
    subplot(3,3,k)
    [width, time]=hist(widthall2(:,k),30);
    bar(time,width)
    title(sprintf('Width: %d',k))
    
    [idx,C]=kmeans(widthall2(:,k),2);
    figure(2)
    subplot(3,3,k)
    hold on
%     loglog(widthall2(idx==1,k),baseline2(idx==1),'bo')
%     loglog(widthall2(idx==2,k),baseline2(idx==2),'ro')
     plot(widthall2(idx==1,k),baseline2(idx==1),'bo')
     plot(widthall2(idx==2,k),baseline2(idx==2),'ro')
    
%     for d=1:2
%         if d==1
%             place=find(guess2<=.5 & guess2>=0);
%             plot(widthall(place,k),baseline(place),'gs')
%         elseif d==2
%             place=find(guess2>.7);
%             plot(widthall(place,k),baseline(place),'b*')
% %         elseif d==3
% %             place=find(guess2>.4 & guess2<=.6);
% %             plot(widthall(place,k),baseline(place),'r.')
% %         elseif d==4
% %             place=find(guess2<0);
% %             plot(widthall(place,k),baseline(place),'rd')
%         end
%    end
%     vectG=linspace(0,max(widthall(:,k)),1000); 
%     plot(vectG,slope(k,1).*vectG+slope(k,2),'k');
    ylabel('Firing Rate (spikes/s)')
    xlabel('Width (ms)')
    %axis([0 max(widthall(:,k)) 0 26])
        title(sprintf('Width: %d',k))
end

figure;
hist(baseline2,30)
xlabel('Firing Rate')

% save(filename,'s')
% save(file_slope,'slope')
%save('\widthclass.mat','spks','Erinn', 'widthall','guess2','baseline','temp','spksvar','axon')