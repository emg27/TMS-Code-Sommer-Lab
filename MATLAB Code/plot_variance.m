%%Plot the average waveform and its variance.
function []=plot_variance(x,lower,upper,color)%,waveplot)
%figure(waveplot)
set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);
end

