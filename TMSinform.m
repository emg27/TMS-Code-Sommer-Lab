load('SFN_TMS.mat')

for n=1:size(s,2)
    name=s(n).Name;
    if size(name,2)>0
    placement=strmatch(name,FileName);
    s(n).BrainArea=BrainArea(placement);
    s(n).Stim=StimSham(placement);
    s(n).Intensity=Intensity(placement);
    s(n).Locate=GridLocation(placement);
    s(n).Depth=RawDepthmm(placement);
    end
end
    
save('\SFN_TMS.mat','s')