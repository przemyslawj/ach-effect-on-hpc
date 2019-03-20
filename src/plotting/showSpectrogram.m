datarootdir = '/mnt/DATA/Clara/baseline/2018-09-06';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

subplots = 1;

meta = ReadMeta(binName, path);

tracking_filepath = get_trackingfilepath(datarootdir, binName);

binNameParts = strsplit(binName, '_g0');
expname = binNameParts{1};

secondOffset = 0;
lengthSeconds = min(str2double(meta.fileTimeSecs) - secondOffset, 60);

nChans = meta.nChans;

animal_code = binName(1:2);

channelList = findSelectedChannels(...
    '/mnt/DATA/Clara/ymaze/selected_electrodes.csv',...
    animal_code);

if isempty(channelList)
    channelList=1:nChans;
end

fs = 600;
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);
if ~isempty(time_mouse_arrived)
    time_mouse_arrived = time_mouse_arrived([2 5 6],:);
end

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
%dataArray = filter50Hz(dataArray, fs);

%% Filter LFP for ripples
filtered = zeros(size(dataArray));
passband = [100 250];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);

for chan = 1:nChans
    % Use MATLAB built-in function
    filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
end

freqrange = 1:2:50;
for i = 1:numel(channelList)
    channel = channelList(i);
    figure('Name', [binName '-channel-' num2str(channel)]);
    % Raw signal
    subplot(4,1,1);
    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints, dataArray(channel,:));
    xstd = std(dataArray(channel,:));
    
    subplot(4,1,2);
    %plot(timepoints, filtered(channel,:));
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(timepoints', filtered(channel,:)', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 4.0 0.01],...
                                 'durations', [10 40 350]);
    ripple_starts = [];
    ripple_ends = [];
    if ~isempty(ripples)                               
        ripple_starts = ripples(:,1);
        ripple_ends = ripples(:,3);
    end
    plotSWR(timepoints, filtered(channel,:), fs, ripple_starts, ripple_ends);
    
    % Spectrogram signal
    [wt, wfreqs]=cwt(dataArray(channel,:), 'morse', fs, 'ExtendSignal', true, ...
        'VoicePerOctave', 30, 'WaveletParameters', [3 120]);
    wt_pow = abs(wt).^2;
    low_freqs = find(wfreqs <= 20);
    high_freqs = find(wfreqs < 200);    
    
    if subplots == 0
        figure;
    end
    subplot(4,1,3);
    draw_cwt(wt_pow(high_freqs,:), timepoints, wfreqs(high_freqs));
    draw_keypoints(time_mouse_arrived, [min(wfreqs(high_freqs)), max(wfreqs(high_freqs))], lengthSeconds, secondOffset);
    
    subplot(4,1,4);
    draw_cwt(wt_pow(low_freqs,:), timepoints, wfreqs(low_freqs));
    draw_keypoints(time_mouse_arrived, [min(wfreqs(low_freqs)), max(wfreqs(low_freqs))], lengthSeconds, secondOffset);
    ax = gca;
    ax.XAxis.Visible = 'on';
    xlabel('Time (sec)');  
    
%     figure('Name', ['Pwelch' '-channel-' num2str(channel)]);
%     [pxx, freqs] = pwelch(dataArray(channel,1:fs*15), ...
%                 floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
%     plot(freqs, 10*log10(pxx))
%     
%     figure('Name', ['CWT' '-channel-' num2str(channel)]);
%     [wt, wfreqs]=cwt(dataArray(channel,fs*1:fs*15), 'amor', fs);
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

function [] = draw_cwt(cfs,time,freq)
    A = z_score(cfs);
    args = {time, freq, A};
    surf(args{:},'edgecolor','none');
    view(0,90);
    axis xy;
    axis tight;
    ax = gca;
    ax.XAxis.Visible = 'off'; % remove x-axis
    %ax.YTick = 0:10:200;
    
    shading interp; 
    %colormap(parula(256));
    colormap(jet);
    colorbar('off');
    %h = colorbar;
    %h.Label.String = 'z-score';
    ylabel('Frequency (Hz)');
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
