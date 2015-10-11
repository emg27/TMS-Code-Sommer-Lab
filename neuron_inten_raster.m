%%Create a raster plot for a single cell example that builds up in
%%intensities. The code will work by first searching the structure of
%%collecting data you are interested in. The code will then arrange a
%%raster plot of the neuron in the order of intensities for some N number
%%of first TMS pulses. Additionally each plot will be a slightly different
%%color in order to make it easier to see the shift in intensity responses.

%Clear and set up the workspace
% clear
% [filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
% load([pathname filename])

%Establish your desire values
date='20130828';%'20150810'; %Date to take the spiking data
stimCond='Stim'; %Set which stimulation condition you are interested in
neuron=9; %Which cluster will the neuron be in
N=10; %The number of TMS pulses (trials you would like to show on the plot)
ta=4000;
tb=4000;

%Function
ColorSet=@(targetcode, shift) [1/targetcode targetcode/9 shift];

%Find your file locations
fileloc=find(arrayfun(@(n) strncmp(s(n).Name, '20150810',8), 1:numel(s)));
if size(fileloc,2)<=0
    break
end

raster=figure;
%line([0 0],[0 9*N],'Color','k')
line([-tb ta],[10 10;20 20;30 30;40 40;50 50;60 60;70 70;80 80;90 90],...
    'Color', 'k', 'LineStyle', '--')
check=zeros(1,9);
for n=1:size(fileloc,2)
    block=s(fileloc(n));
    inten=str2num(block.Intensity{1})/10;
    if strcmp(block.Stim(1),stimCond) & length(s(n).Pulses)>0 & check(inten)==0
        clusterpos=find(block.clusters==neuron);
    else
        continue
    end
    Pulses=s(n).Pulses(1:N);
    
    shift=(inten-1)*N;
    figure(raster)
    color=Raster(Pulses,tb,ta,1000*block.times(clusterpos),shift);
    
    set(color,'Color',ColorSet((9-inten)+1,.5));
    check(inten)=1;
end