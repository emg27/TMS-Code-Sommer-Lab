%%plot the baseline before pulse through time during a tms recording 

%Clear workspace
close all

%load data
file='base_save_SFN.mat';
[FileName,PathName] = uiputfile('*.mat');
load(file);
load([PathName,FileName])

%set the parameters
t_period=500; %ms
base_bef=[];
base_aft=[];

dose=figure;

for n=1:size(base_save(1,:))
end