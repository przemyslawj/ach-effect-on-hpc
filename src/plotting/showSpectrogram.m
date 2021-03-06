datarootdir = '/mnt/DATA/Clara/baseline/2018-09-06';
datarootdir = '/mnt/DATA/chat_ripples/y-maze/2019-11-13';
%datarootdir = '/mnt/DATA/chat_ripples/y-maze/2020-08-28';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

subplots = 1;
selected_channels_only = 1;
use_diode = 1;
reversed_channel_map = 1;

meta = ReadMeta(binName, path);

tracking_filepath = get_trackingfilepath(datarootdir, binName);

binNameParts = strsplit(binName, '_g0');
expname = binNameParts{1};

secondOffset = 2.5;
lengthSeconds = min(str2double(meta.fileTimeSecs) - secondOffset, 50);

animal_code = binName(1:2);
channelTable = readChannelTable(...
    '/mnt/DATA/chat_ripples/channel_desc/channels.csv',...
    animal_code, meta, reversed_channel_map, selected_channels_only, use_diode);
fs = 600;
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);
if ~isempty(time_mouse_arrived)
    time_mouse_arrived = time_mouse_arrived([2 5 6],:);
end

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);
dataArray = dataArray(channelTable.rec_order,:);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
%dataArray = filter50Hz(dataArray, fs);

if use_diode
    [dataArray, channelTable] = subtractDiodeSignal(dataArray, channelTable);
end
%% Filter LFP for ripples
filtered = applyRippleFilter(dataArray, channelTable, fs);

%lengthSeconds = 0.5;
%startSec = 28.7;
%timeIndecies = (startSec * fs) : ((startSec + lengthSeconds) * fs);
timeIndecies = 1:size(dataArray, 2);


for chan_i = 1:size(channelTable, 1)
    loc = channelTable.location(chan_i);
    if strcmp(loc, 'EMG') || strcmp(loc, 'Laser')
        continue
    end
    figure('Name', [binName '-' channelTable.channel_name{chan_i} ...
                    ', channel:' num2str(channelTable.channel(chan_i))]);
    % Workaround to force saving as svg with all the elements
    %set(0, 'DefaultFigureRenderer', 'painters');
    
    % Raw signal
    subplot(5,1,1);
    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints(timeIndecies), dataArray(chan_i,timeIndecies));
    xstd = std(dataArray(chan_i,:));
    
    subplot(5,1,2);
    %plot(timepoints, filtered(channel,:));
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(timepoints',...
                                 filtered(chan_i,:)', ...
                                 'frequency', fs, ...
                                 'thresholds', [2.5 6.0 0.01],...
                                 'durations', [10 40 350]);
    ripple_starts = [];
    ripple_ends = [];
    if ~isempty(ripples)                               
        ripple_starts = ripples(:,1);
        ripple_ends = ripples(:,3);
    end
    plotSWR(timepoints(timeIndecies), filtered(chan_i,timeIndecies), fs, ripple_starts, ripple_ends);
    
    % Spectrogram signal
    [wt, wfreqs]=cwt(dataArray(chan_i,:), 'morse', fs, 'ExtendSignal', true, ...
        'VoicePerOctave', 30, 'WaveletParameters', [3 120]);
    wt_pow = abs(wt).^2;
    low_freqs = find(wfreqs <= 15);
    high_freqs = find(wfreqs < 250 & wfreqs > 15);
    
    if subplots == 0
        figure;
    end

%     % Workaround to force saving as svg with all the elements
%     set(0, 'DefaultFigureRenderer', 'painters');
    wt_pow_scaled = matScale(wt_pow, 1000);
    scaledTimepoints = matScale(timepoints,1000);
    hAx = [];
    hAx(1) = subplot(5,1,[3 4]);
    A = z_score(wt_pow_scaled(high_freqs,:));
    draw_cwt(A, scaledTimepoints, wfreqs(high_freqs));
    h = colorbar;
    zlimits = h.Limits;
    %draw_keypoints(time_mouse_arrived, [min(wfreqs(high_freqs)), max(wfreqs(high_freqs))], lengthSeconds);
    
    hAx(2) = subplot(5,1,5);
    A = zscore(wt_pow_scaled(low_freqs,:));
    draw_cwt(A, scaledTimepoints, wfreqs(low_freqs));
    h = colorbar;
    h.Label.String = 'z-score';
    h.Limits = zlimits;
    %draw_keypoints(time_mouse_arrived, [min(wfreqs(low_freqs)), max(wfreqs(low_freqs))], lengthSeconds);
    ax = gca;
    ax.XAxis.Visible = 'on';
    xlabel('Time (sec)');
    linkaxes(hAx,'x');
    
%     figure('Name', ['Pwelch' '-channel-' num2str(chan_i)]);
%     [pxx, freqs] = pwelch(dataArray(chan_i,1:fs*15), ...
%                 floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
%     plot(freqs, 10*log10(pxx))
%     
%     figure('Name', ['CWT' '-channel-' num2str(chan_i)]);
%     [wt, wfreqs]=cwt(dataArray(chan_i,fs*1:fs*15), 'amor', fs);
%     wt_pow = abs(wt).^2;
%     plot(wfreqs, median(wt_pow, 2))
end

function A = z_score(cfs)
    A = cfs - mean(cfs, 2);
    A = bsxfun(@rdivide, A, std(A, [], 2));
    maxZscore = ones(size(A)) * 6;
    minZscore = ones(size(A)) * -3;
    A = bsxfun(@min, A, maxZscore);
    A = bsxfun(@max, A, minZscore);
end


function [] = draw_keypoints(time_mouse_arrived, ylim, lengthSeconds)
    if isempty(time_mouse_arrived)
        return
    end
    for i = 1:numel(time_mouse_arrived.sec)
        x = time_mouse_arrived.sec(i);
        if x < lengthSeconds && x > 0
            line([x, x], ylim, 'Color', 'white');
        end
    end
end

%% Utilify functions
function [scaledM] = matScale(M, required_width)
    % Stretch the trial period length
    mwidth = size(M, 2);
    if mwidth < required_width
        M = repelem(M, 1, ceil(required_width / mwidth));
    end

    % Compress the trial period
    blockSize = [1, floor(mwidth / required_width)];
    meanFilterFun = @(blockStruct) mean2(blockStruct.data(:));
    scaledM = blockproc(M, blockSize, meanFilterFun);
    scaledM = scaledM(:, 1:required_width);
end
