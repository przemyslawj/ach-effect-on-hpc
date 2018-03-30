path = '/media/prez/DATA/Prez/N&A_rest/';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

secondOffset = 0;
meta = ReadMeta(binName, path);

%lengthSeconds = min(40, str2double(meta.fileTimeSecs) - secondOffset);
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
%lengthSeconds = 3;
fs = meta.nSamp;
%nChans = str2double(meta.nChans);
nChans=32;
channelList = 1:nChans;
channelList = [1 8];

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 1250;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

dataArray = filter50Hz(dataArray, fs);

%% Filter LFP
filtered = zeros(size(dataArray));
passband = [100 250];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);

for chan = 1:nChans
    %filtered(i,:) = FiltFiltM(b, a, dataArray(i,:));
    % Use MATLAB built-in function
    filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
end

%% Spike detection
spikeFiltered = zeros(size(dataArray));
highpass = 300;
filterOrder = 9;
filterRipple = 40;
[b, a] = cheby2(filterOrder, filterRipple, highpass / nyquist, 'high');

for chan = 1:nChans
    %spikeFiltered(i,:) = FiltFiltM(b, a, dataArray(i,:));
    spikeFiltered(chan,:) = filtfilt(b, a, dataArray(chan,:));
end
spikeTimes = spike_amp_detect(dataArray(chan,:), spikeFiltered(chan,:));

%% Detect spikes

% Remove population spikes as cable movement artifacts
% (>=25% channels firing during the same time bin) 
binsize = round(0.004 * fs);
binranges = 1 : binsize : (size(spikeFiltered, 2) + binsize);
binnedspikes = zeros(nChans, length(binranges));

spikes_i =  zeros(size(dataArray));
spike_thresholds = zeros(1,nChans);
for chan = 1:nChans
   [~, spike_indecies, spike_thr] = spike_amp_detect(dataArray(chan,:), spikeFiltered(chan, :));
   spike_thresholds(chan) = spike_thr;
   spikes_i(chan, spike_indecies) = 1;
   binnedspikes(chan,:) = histc(spike_indecies, binranges);
end

chansFiring = sum(binnedspikes) / nChans;
keptSpikeBins = chansFiring <= 0.25;
for start_time = 1 : binsize : (size(spikes_i,2) - binsize)
    times = start_time:(start_time+binsize-1);
    bin_i = ceil(start_time / binsize);
    spikes_i(:,times) = spikes_i(:,times) * keptSpikeBins(bin_i);
end

%% Plot SWR and spikes
for i = 1:numel(channelList)
    chan = channelList(i);
    time=(1:size(filtered,2))/fs;
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(chan,:)', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 4 0.02],...
                                 'durations', [30 20 300]);
    spikeTimes = find(spikes_i(chan,:) > 0) / fs;
    if ~isempty(ripples)                               
        ripple_starts = ripples(:,1);
        ripple_ends = ripples(:,3);
    else
        ripple_starts = [];
        ripple_ends = [];
    end 
    fprintf('Channel %d: %d ripples, %d spikes\n', chan, ...
        length(ripple_starts), length(spikeTimes));
    
    figure('name', sprintf('Channel %d', chan));
    plotData(fs, dataArray(chan,:), filtered(chan,:), spikeFiltered(chan,:), ...
        ripple_starts, ripple_ends, normalizedSquaredSignal', ...
        spikeTimes, spike_thresholds(chan), 0);
    
    sliderVars = struct('fs', fs, 'lengthSeconds', lengthSeconds,...
        'data', dataArray(chan,:), 'filtered', filtered(chan,:), ...
        'spikeFiltered', spikeFiltered(chan,:), ...
        'spikeTimes', spikeTimes, ...
        'spikeThreshold', spike_thresholds(chan), ...
        'ripple_starts', ripple_starts, 'ripple_ends', ripple_ends,...
        'normalizedSquaredSignal', normalizedSquaredSignal');
    sliderHandle = uicontrol('Style', 'slider', ...
              'Position', [10 20 500 20]); 
    set(sliderHandle,'Callback',{@sliderCallback, sliderVars});
    
end


