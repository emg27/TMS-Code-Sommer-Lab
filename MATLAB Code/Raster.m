function[points,position, num_spike,blocks]=Raster(pulsetime,time_b,time_a,spikes,blocks,num_spike)
%%Raster Plots: Created March 14, 2013 This file will create basic raster
%%plots for the specified for the values necessary:
%%Raster(pulsetime,time_b,time_a,spikes,blocks,num_spikes)
%%Pulsetime= The time that the tms pulse begin, this will be the center 
%%marker of the raster plot 
%%time_b= The amount of time before the pulse to include in the raster
%%time_a= The amount of time after the pulse to include in the raster
%%spikes= The spike times, this should be a 2 column matrix, where the 
%%second column contains the spike times and the first has the cluster
%%number
%%blocks= The number blocks per each trial

%%function
if nargin==4
    blocks=zeros(length(pulsetime),1);
    num_spike=zeros(length(pulsetime),1);
end
points=[];
position=[];
for n=1:length(pulsetime)
    begin=pulsetime(n)-time_b; %This is the time marker for the min time before each pulse
    fin=pulsetime(n)+time_a; %This is the time marker for the min time before each pulse
    placemark=find(spikes>=begin & spikes<=fin); %Finds the times positions in between the pulses
    trial=spikes(placemark); %Saves the exact time stamp
    num_spike(n)=num_spike(n)+length(placemark);
    if(num_spike(n)>0)
        blocks(n)=blocks(n)+1;
    end
    
   point=trial-pulsetime(n);
   if ~isempty(point)
       pointA=line([point point]', [n-0.9 n-0.1]); hold on
       points=[points pointA'];
   end
    %num_plot(n,1)=text(-time_b-20,[n-.5],sprintf('n=%0.0d',num_spike(n))); %prints number of spikes per trial
    %if(max(blocks)>1)
    %    num_plot(n,2)=text(time_a+5,[n-.5],sprintf('b=%d',blocks(n))); %prints number of blocks per trial
    %end
  
position=[position; placemark];
end

    end