%%Plot the entire raster over an entire block and the waveformdlkflkjfdl
%%blockraster(alignmenttimes,times,placement,color)
%%alignments: The time when the tms pulse happened
%%times: The times when the spikes occurred
%%placement: Where you want the block plot to be placed on a figure, ie if
%%you want the block to be placed on y=7 or y=-2.
%%color: What color you want the spike rasters to be
function[]=blockraster(alignmenttimes,times,placement,color)%,spikes,blockplot,waveplot)


plot([alignmenttimes alignmenttimes],[-1.5+placement 1.5+placement],'k')
hold on;
plot([times times]',[-.5+placement .5+placement],'Color',color)

end
