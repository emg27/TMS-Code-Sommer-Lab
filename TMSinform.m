%load('SFN_TMS.mat')
% [File,PathName] = uiputfile('*.mat')
% load([PathName File])

%%Create FileName 
dash=strfind(Date,'/');

for k=1:size(Date,1)
    dashpos=dash{k};
    date=Date{k};
    block=BlockNumber{k};
    if isnumeric(block)
        block=num2str(block);
    end
    if size(dashpos,2)<2
        FileName(k)=cellstr('');
    elseif dashpos(1)==2 
        if diff(dashpos)==3
            FileName(k)=cellstr([date(end-3:end) '0' date([1 3 4]) '_' block]);
        elseif diff(dashpos)==2
            FileName(k)=cellstr([date(end-3:end) '0' date(1) '0' date(3) '_' block]);
        end
    elseif dashpos(1)==3
        if diff(dashpos)==3
            FileName(k)=cellstr([date(end-3:end) date([1 2 4 5]) '_' block]);
        elseif diff(dashpos)==2
            FileName(k)=cellstr([date(end-3:end) date([1 2]) '0' date(4) '_' block]);
        end
    end
end
FileName=FileName';

for n=1:size(s,2)
    name=s(n).Name;
    if size(name,2)>0
    placement=strmatch(name,FileName);
    s(n).BrainArea=BrainArea(placement);
    s(n).Stim=StimSham(placement);
    s(n).Intensity=Intensity(placement);
    s(n).Locate=GridLocation(placement);
    %s(n).Depth=RawDepthmm(placement);
    end
end

[File,PathName] = uiputfile('*.mat','Save Data File')
save([PathName File],'s')