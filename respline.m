function [w1] = respline(wav, freq, pts, tmin, tmax)
% Function takes a waveform, creates an interpolated cubic spline, then
% reextracts that spline to a new set of data points at a desired frequency
%% Define Parameters
% Time vector for arbitrary wave
twav = linspace(0, (length(wav)-1)/freq, length(wav));
% Time vector for desired wave
tt = linspace(tmin, tmax, pts);
% time between each sample of desired wave
delT = (abs(tmin) + abs(tmax))/pts;
% time vector for realignment
ttAlign = linspace(0, max(twav), ceil(max(twav)/delT));

%% Take Alignment Spike
wwAlign = spline(twav, wav, ttAlign);


%% Exceptions


%% 


end