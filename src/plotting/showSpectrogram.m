%datarootdir = '/mnt/DATA/Clara/baseline/2018-09-06';
datarootdir = '/mnt/DATA/chat_ripples/y-maze/2019-11-13';
%datarootdir = '/mnt/DATA/chat_ripples/sleep/2019-12-03';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

subplots = 1;
selected_channels_only = 1;
use_diode = 0;

meta = ReadMeta(binName, path);

tracking_filepath = get_trackingfilepath(datarootdir, binName);

binNameParts = strsplit(binName, '_g0');
expname = binNameParts{1};

secondOffset = 3;
lengthSeconds = min(str2double(meta.fileTimeSecs) - secondOffset, 40);

animal_code = binName(1:2);
channelTable = readChannelTable(...
    '/mnt/DATA/chat_ripples/channel_desc/channels_reversed.csv',...
    animal_code, meta, selected_channels_only, use_diode);
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
    % Raw signal
    subplot(4,1,1);
    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints(timeIndecies), dataArray(chan_i,timeIndecies));
    xstd = std(dataArray(chan_i,:));
    
    subplot(4,1,2);
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
    low_freqs = find(wfreqs <= 20);
    high_freqs = find(wfreqs < 200 & wfreqs > 20);
    
    if subplots == 0
        figure;
    end
    subplot(4,1,3);
    A = z_score(wt_pow(high_freqs,timeIndecies));
    draw_cwt(A, timepoints(timeIndecies), wfreqs(high_freqs));
    %draw_keypoints(time_mouse_arrived, [min(wfreqs(high_freqs)), max(wfreqs(high_freqs))], lengthSeconds, secondOffset);
    
    subplot(4,1,4);
    A = zscore(wt_pow(low_freqs,timeIndecies));
    draw_cwt(A, timepoints(timeIndecies), wfreqs(low_freqs));
    %draw_keypoints(time_mouse_arrived, [min(wfreqs(low_freqs)), max(wfreqs(low_freqs))], lengthSeconds, secondOffset);
    ax = gca;
    ax.XAxis.Visible = 'on';
    xlabel('Time (sec)');  
    
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


function [] = draw_keypoints(time_mouse_arrived, ylim, lengthSeconds, secondOffset)
    if isempty(time_mouse_arrived)
        return
    end
    for i = 1:numel(time_mouse_arrived.sec)
        x = time_mouse_arrived.sec(i) - secondOffset;
        if x < lengthSeconds && x > 0
            line([x, x], ylim, 'Color', 'white');
        end
    end
end
