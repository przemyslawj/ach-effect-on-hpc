path = '/mnt/DATA/chat_ripples/baseline/';
%path = '/mnt/DATA/Prez/N&A_rest/2018-03-01/signal';
%path = '/mnt/DATA/chat_ripples/y-maze/2019-11-13/signal';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

selected_channels_only = 0;
use_diode = 1;
reversed_channel_map = 1;

secondOffset = 0;
meta = ReadMeta(binName, path);
lengthSeconds = min(180, str2double(meta.fileTimeSecs) - secondOffset);
%lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
%lengthSeconds = 30;

uppercase_idx = find(isstrprop(binName, 'upper'));
animal_code = binName(uppercase_idx(1:2));
%'/mnt/DATA/chat_ripples/channel_desc/channels_reversed_baseline.csv',...
channelTable = readChannelTable(...
    '/mnt/DATA/chat_ripples/channel_desc/channels.csv',...
    animal_code, meta, reversed_channel_map, selected_channels_only, use_diode);

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 1250;
%dataArray = filter50Hz(dataArray, fs);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
dataArray = dataArray(channelTable.rec_order,:);
time = (1:size(dataArray,2)) / fs;
if use_diode
    [dataArray, channelTable] = subtractDiodeSignal(dataArray, channelTable);
end
filtered = applyRippleFilter(dataArray, channelTable, fs);
emgIdx = find(strcmp(channelTable.location, 'EMG'));

keep_sample_fewer = excludeEMGNoisePeriods(dataArray(emgIdx,:), fs * 1);
%% Plot SWR
for chan = 1:size(dataArray, 1)
    ripple_detection_signal = GetRippleSignal(filtered(chan, :)', fs);
    std_estimate = std(ripple_detection_signal(keep_sample_fewer));
    fprintf('filtered std %.8f\n', std(filtered(chan, keep_sample_fewer)));
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(chan,:)', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 6 0.01],...
                                 'durations', [10 20 300],...
                                 'std', std_estimate);
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
    rec_len_sec = 4;
    plotData(fs, dataArray(chan,:), filtered(chan,:), ...
        ripple_starts, ripple_ends, normalizedSquaredSignal', ...
        0, rec_len_sec);

    sliderVars = struct('fs', fs, 'lengthSeconds', lengthSeconds,...
        'data', dataArray(chan,:), 'filtered', filtered(chan,:), ...
        'ripple_starts', ripple_starts, 'ripple_ends', ripple_ends,...
        'normalizedSquaredSignal', normalizedSquaredSignal',...
        'slider_window_len_sec', 2 * rec_len_sec);
    sliderHandle = uicontrol('Style', 'slider', ...
              'Position', [10 20 1000 20], ...
              'UserData', struct('rec_len_sec', rec_len_sec, 'slider_pos', 0),...
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


