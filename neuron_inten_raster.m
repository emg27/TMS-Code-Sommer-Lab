%%Create a raster plot for a single cell example that builds up in
%%intensities. The code will work by first searching the structure of
%%collecting data you are interested in. The code will then arrange a
%%raster plot of the neuron in the order of intensities for some N number
%%of first TMS pulses. Additionally each plot will be a slightly different
%%color in order to make it easier to see the shift in intensity responses.

% %Clear and set up the workspace
% clear
% [filename, pathname]=uigetfile('*.mat')%'oxford_2014.mat';
% load([pathname filename])
% close all

%Establish your desire values
date='20150109';%'20141009';
notall=0;
%'20150123';

%'20150123';%'20140625';%'20150810'; %Date to take the spiking data
stimCond='Stim'; %Set which stimulation condition you are interested in
neuron=1; %Which cluster will the neuron be in
N=10; %The number of TMS pulses (trials you would like to show on the plot)
ta=300;
tb=100;

%Function
ColorSet=@(targetcode, shift) [1/targetcode targetcode/9 shift];

%Find your file locations
fileloc=find(arrayfun(@(n) strncmp(s(n).Name, date,8), 1:numel(s)));
if size(fileloc,2)<=0
    break
end

wave=figure;
raster=figure;

line([0 0],[0 9*N],'Color','k')
line([-tb ta],[N N;2*N 2*N;3*N 3*N;4*N 4*N;5*N 5*N;6*N 6*N;7*N 7*N;8*N 8*N;9*N 9*N],...
    'Color', 'k', 'LineStyle', '--')
check=zeros(1,9);
k=1;
waves=[]
for n=1:size(fileloc,2)
    block=s(fileloc(n));
    inten=str2num(block.Intensity{1})/10;
    if strcmp(block.Stim(1),stimCond) & length(block.Pulses)>0 & inten<=9 & check(inten)==0
        clusterpos=find(block.clusters==neuron);
    else
        continue
    end
    checkpos=0;
    if notall==1
        if inten>=6
            checkpos=2;
        elseif inten>=4
            checkpos=1;
        end
    end
    Pulses=block.Pulses(1:N);
    
    shift=(inten-(1+checkpos))*N;
    figure(raster)
    color=Raster(Pulses,tb,ta,1000*block.times(clusterpos),shift);
    
    set(color,'Color',ColorSet((9-inten)+1,.5));
    check(inten)=1;
    
    figure(wave)
    subplot(3,3,k)
    meanWor=mean(block.waveforms(clusterpos,:));
    stdWor=std(block.waveforms(clusterpos,:));
    t1= 1000*(0:1:length(meanWor)-1)/block.FireRate;
    time=linspace(min(t1),max(t1),1000);
    meanW=spline(t1,meanWor,time);
    stdW=spline(t1,stdWor,time);
    plot_variance(time,meanW+stdW,meanW-stdW,'b');
    hold on;
    plot(time,meanW,'k');
    axis([0 time(end) -.15 .15])
    title(num2str(inten*10))
    
    k=k+1;
    %waves=[waves; block.waveforms(clusterpos,:)];
end
figure(raster)
title('Monkey M Single Cell Example')
ylabel('Intensity Trials')
xlabel('Time (ms)')
if notall==1
    ylim([0 70])
elseif notall==2
    ylim([0 80])
end