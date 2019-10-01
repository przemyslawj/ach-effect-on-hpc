path = '/mnt/DATA/Clara/ymaze/2018-08-17/signal';
path = '/mnt/DATA/Prez/N&A_rest/2018-03-01/signal';
path = '/mnt/DATA/Clara/diode_baseline/20190906/BD031ActiveStimON1_g0';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

secondOffset = 0;
meta = ReadMeta(binName, path);
%lengthSeconds = min(40, str2double(meta.fileTimeSecs) - secondOffset);
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
%lengthSeconds = 30;
animal_code = binName(1:2);
channelTable = findOmneticChannels(...
    '/mnt/DATA/Clara/diode_baseline/channels.csv',...
    animal_code);
keepChannels = intersect(meta.snsSaveChanSubset, channelTable.channel, 'sorted');
channelTable = channelTable(ismember(channelTable.channel, keepChannels), :);

nChans = size(channelTable, 1);

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = meta.nSamp / 10;
%dataArray = filter50Hz(dataArray, fs);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
dataArray = dataArray(ismember(meta.snsSaveChanSubset, keepChannels),:);

%% Filter LFP
filtered = zeros(size(dataArray));
passband = [80 250];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);

for chan = 1:nChans
    filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
end

%% Plot SWR
for chan = 1:nChans
    time=(1:size(filtered,2))/fs;
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(chan,:)', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 6 0.01],...
                                 'durations', [10 40 350]);
    if ~isempty(ripples)                               
        ripple_starts = ripples(:,1);
        ripple_ends = ripples(:,3);
    else
        ripple_starts = [];
        ripple_ends = [];
    end 
    fprintf('Channel %s: %d ripples\n', channelTable.channel_name{chan}, ...
        length(ripple_starts));
    
    figure('name', sprintf('Channel %s', channelTable.channel_name{chan}),...
           'Position', [100 1000 1100 600]);
    nfigure = length(findobj('type','figure'));
    plotData(fs, dataArray(chan,:), filtered(chan,:), ...
        ripple_starts, ripple_ends, normalizedSquaredSignal', ...
        0, 2);
    
    sliderVars = struct('fs', fs, 'lengthSeconds', lengthSeconds,...
        'data', dataArray(chan,:), 'filtered', filtered(chan,:), ...
        'ripple_starts', ripple_starts, 'ripple_ends', ripple_ends,...
        'normalizedSquaredSignal', normalizedSquaredSignal');
    sliderHandle = uicontrol('Style', 'slider', ...
              'Position', [10 20 1000 20], ...
              'UserData', struct('rec_len_sec', 2, 'slider_pos', 0),...
              'SliderStep', [1 / lengthSeconds / 4 , 10 / lengthSeconds / 4],...
              'Tag', sprintf('timeSlider_%d', nfigure));
    set(sliderHandle,'Callback',{@sliderCallback, sliderVars});
    
    timewindowHandle = uicontrol('Style', 'slider', ...
          'Position', [10 50 100 20],...
          'SliderStep', [0.1, 0.3],...
          'Tag', sprintf('timewindowSlider_%d', nfigure));
    timewindowHandle.Value = 0.5;
    set(timewindowHandle,'Callback',{@timewindowCallback, sliderVars});
end


