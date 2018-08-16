function [ripples,sd,normalizedSquaredSignal] = MyFindRipples(time, signal, varargin)

%FindRipples - Find hippocampal ripples.
%
%  USAGE
%
%    [ripples,stdev,noise] = FindRipples(filtered,<options>)
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.
%
%    filtered       ripple-band filtered LFP (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5])
%     'durations'   minimum inter-ripple interval, and minimum and maximum
%                   ripple durations, in ms (default = [30 20 100])
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'show'        plot results (default = 'off')
%     'noise'       noisy ripple-band filtered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%    =========================================================================
%
%  OUTPUT
%
%    ripples        for each ripple, [start_t peak_t end_t peakNormalizedPower]
%    stdev          standard deviation of the NSS (can be reused subsequently)
%    noise          ripple-like activity recorded simultaneously on the noise
%                   channel (for debugging info)
%
%  SEE
%
%    See also FilterLFP, RippleStats, SaveRippleEvents, PlotRippleStats.

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
frequency = 1250;
show = 'off';
lowThresholdFactor = 2; % Ripple envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 5; % Ripple peak must exceed highThresholdFactor*stdev
minInterRippleInterval = 30; % in ms
minRippleDuration = 30; % in ms
maxRippleDuration = 250; % in ms
lowThresholdMV = 0.02;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin)
	if ~ischar(varargin{i})
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
	end
	switch(lower(varargin{i}))
		case 'thresholds'
			thresholds = varargin{i+1};
% 			if ~isdvector(thresholds,'#2','>0'),
% 				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
% 			end
			lowThresholdFactor = thresholds(1);
			highThresholdFactor = thresholds(2);
            lowThresholdMV = thresholds(3);
		case 'durations'
			durations = varargin{i+1};
			if length(durations) == 2
				minInterRippleInterval = durations(1);
				maxRippleDuration = durations(2);
			else
				minInterRippleInterval = durations(1);
				minRippleDuration = durations(2);
				maxRippleDuration = durations(3);
			end
		case 'frequency'
			frequency = varargin{i+1};
% 			if ~isdscalar(frequency,'>0'),
% 				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
% 			end
		case 'show'
 			show = varargin{i+1};
% 			if ~isstring(show,'on','off'),
% 				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
% 			end
		case 'noise'
			noise = varargin{i+1};
			if ~isdmatrix(noise) | size(noise,1) ~= size(filtered,1) | size(noise,2) ~= 2
				error('Incorrect value for property ''noise'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        otherwise
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
	end
end

% Parameters
windowLength = round(frequency/11);
% make window of odd length
windowLength = floor(windowLength / 2) * 2 + 1; 

squaredSignal = signal.^2;
window = ones(windowLength,1)/windowLength;

[normalizedSquaredSignal,sd] = unity(Filter0(window, squaredSignal));

% Detect ripple periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1
	start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if ~isempty(stop) && ~isempty(start) && start(1) > stop(1)
	stop(1) = [];
	start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass)
	disp('No ripples detected');
    ripples = [];
    bad = [];
    sd = [];
	return
else
	%disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge ripples if inter-ripple period is too short
minInterRippleSamples = minInterRippleInterval/1000*frequency;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - ripple(2) < minInterRippleSamples
		% Merge
		ripple = [ripple(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; ripple];
		ripple = firstPass(i,:);
	end
end
secondPass = [secondPass ; ripple];
if isempty(secondPass)
	disp('Ripple merge failed');
    ripples = [];
    bad = [];
    sd = [];
	return
else
	%disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
    maxAbsVal = max(abs(signal([secondPass(i,1):secondPass(i,2)])));
	if maxValue > highThresholdFactor && maxAbsVal > lowThresholdMV
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end
if isempty(thirdPass)
	%disp('Peak thresholding failed.');
    ripples = [];
    bad = [];
    sd = [];
	return
else
	%disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end

% Detect negative peak position for each ripple
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1)
	[minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end

% Discard ripples that are too short or too long
ripples = [time(thirdPass(:,1)) time(peakPosition) time(thirdPass(:,2)) peakNormalizedPower];
duration = ripples(:,3)-ripples(:,1);
ripples(duration < minRippleDuration/1000,:) = [];
ripples(duration > maxRippleDuration/1000,:) = [];
%disp(['After min and max duration test: ' num2str(size(ripples,1)) ' events.']);


% Optionally, plot results
if strcmp(show,'on')
	figure;
    MultiPlotXY([time signal],[time squaredSignal],[time normalizedSquaredSignal]);
    nPlots = 3;
    subplot(nPlots,1,3);
    ylim([0 highThresholdFactor*1.1]);
	for i = 1:nPlots
		subplot(nPlots,1,i);
		hold on;
  		yLim = ylim;
		for j=1:size(ripples,1)
			plot([ripples(j,1) ripples(j,1)],yLim,'g-');
			plot([ripples(j,2) ripples(j,2)],yLim,'k-');
			plot([ripples(j,3) ripples(j,3)],yLim,'r-');
			if i == 3
				plot([ripples(j,1) ripples(j,3)],[ripples(j,4) ripples(j,4)],'k-');
			end
        end
		if mod(i,3) == 0
			plot(xlim,[lowThresholdFactor lowThresholdFactor],'k','linestyle','--');
			plot(xlim,[highThresholdFactor highThresholdFactor],'k-');
		end
    end
end
%disp(['Found: ' num2str(length(secondPass)) ' ripple events.']);

function y = Filter0(b,x)
% Moving average filter calculated with kernel b
    if mod(length(b),2) ~= 1
        error('kernel length should be odd');
    end

    shift = (length(b)-1)/2;

    [y0 z] = filter(b,1,x);

    y = [y0(shift+1:end,:) ; z(1:shift,:)];
end

function [U,stdA] = unity(A)
	meanA = mean(A);
	stdA = std(A);
    U = (A - meanA)/stdA;
end


end

