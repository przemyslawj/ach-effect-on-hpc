path = '/mnt/DATA/Clara/urethane/30s_light/';
[binName, path] = uigetfile('*.txt', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

secondOffset = 0;

fs = 1000;
nChans = 1;
channelList = 1;

datatable = readtable([path filesep binName]);
dataArray = datatable{:,1};
dataArray = dataArray';
lengthSeconds = size(dataArray, 2) / fs;

%dataArray = filter50Hz(dataArray, fs);

%% Filter LFP
filtered = zeros(size(dataArray));
passband = [80 250];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);

for chan = 1:nChans
    %filtered(i,:) = FiltFiltM(b, a, dataArray(i,:));
    % Use MATLAB built-in function
    filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
end

spikeFiltered = zeros(size(dataArray));
spike_thresholds = zeros(1,nChans);

%% Plot SWR and spikes
for i = 1:numel(channelList)
    chan = channelList(i);
    time=(1:size(filtered,2))/fs;
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(chan,:)', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 3 0.01],...
                                 'durations', [10 40 350]);
    spikeTimes = [];
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
    nfigure = length(findobj('type','figure'));
    plotData(fs, dataArray(chan,:), filtered(chan,:), spikeFiltered(chan,:), ...
        ripple_starts, ripple_ends, normalizedSquaredSignal', ...
        spikeTimes, spike_thresholds(chan), 0, 2);
    
    sliderVars = struct('fs', fs, 'lengthSeconds', lengthSeconds,...
        'data', dataArray(chan,:), 'filtered', filtered(chan,:), ...
        'spikeFiltered', spikeFiltered(chan,:), ...
        'spikeTimes', spikeTimes, ...
        'spikeThreshold', spike_thresholds(chan), ...
        'ripple_starts', ripple_starts, 'ripple_ends', ripple_ends,...
        'normalizedSquaredSignal', normalizedSquaredSignal');
    sliderHandle = uicontrol('Style', 'slider', ...
              'Position', [10 20 500 20], ...
              'UserData', struct('rec_len_sec', 2, 'slider_pos', 0),...
              'Tag', sprintf('timeSlider_%d', nfigure));
    set(sliderHandle,'Callback',{@sliderCallback, sliderVars});
    
    timewindowHandle = uicontrol('Style', 'slider', ...
          'Position', [10 50 100 20],...
          'Tag', sprintf('timewindowSlider_%d', nfigure));
    timewindowHandle.Value = 0.5;
    set(timewindowHandle,'Callback',{@timewindowCallback, sliderVars});
end


