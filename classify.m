function[width,p]=classify(meanw,firerate)
if min(meanw)>-.3
%%Waveform classifications

time=1000*(0:1:length(meanw)-1)/firerate;

%%Find the positions of the points
a2=find(meanw==min(meanw));
p2=a2(1);
a4=find(meanw(p2:end)==max(meanw(p2:end)))+p2-1;
p4=a4(1);
a3=find(abs(meanw(p2:p4))==min(abs(meanw(p2:p4))))+p2-1;
p3=a3(1);
a6=find(meanw(p4:end)==min(meanw(p4:end)))+p4-1;
p6=a6(1);
a5=find(abs(meanw(p4:p6))==min(abs(meanw(p4:p6))))+p4-1;
p5=a5(1);

%%Find p1
pointer=.1*meanw(p2);
minref=min(meanw(find(meanw(1:p2)>pointer)));
maxref=max(meanw(find(meanw(1:p2)<=pointer)));

if length(maxref)<=0 | (minref-pointer) < (maxref-pointer)
    p1=find(meanw(1:p2)==minref)-1;
    p1=p1(end);
else
    p1=find(meanw(1:p2)==maxref)-1;
    p1=p1(end);
end

if p1==0
    p1=1;
end
%%Find p7
pointer=.1*meanw(p6);
minref=min(meanw(find(meanw(p6:end)>pointer)+p6-1));
maxref=max(meanw(find(meanw(p6:end)<=pointer)+p6-1));

if length(maxref)<=0 | (minref-pointer) < (maxref-pointer)
    p7=find(meanw(p6:end)==minref)+p6-1;
else
    p7=find(meanw(p6:end)==maxref)+p6-1;
end

%%Width durations
width(1)=time(p4)-time(p1);
width(2)=time(p5)-time(p1);
width(3)=time(p3)-time(p1);
width(4)=time(p4)-time(p2);
width(5)=time(p6)-time(p4);
width(6)=time(p5)-time(p3);
if size(p7)>0 
width(7)=time(p7)-time(p1);
width(8)=time(p7)-time(p5);
else 
    width(7)=NaN
    width(8)=NaN
end
width(9)=time(p6)-time(p2)

p=[p1 p2 p3 p4 p5 p6 p7];
else
    p=NaN(7,1);
    width=NaN(9,1);
end
