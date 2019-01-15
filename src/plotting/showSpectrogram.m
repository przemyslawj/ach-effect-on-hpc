datarootdir = '/mnt/DATA/Clara/ymaze/2018-08-16';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

meta = ReadMeta(binName, path);

tracking_filepath = get_trackingfilepath(datarootdir, binName);

binNameParts = strsplit(binName, '_g0');
expname = binNameParts{1};


secondOffset = 0;
lengthSeconds = min(str2double(meta.fileTimeSecs) - secondOffset, 60);


nChans = meta.nChans;

animal_code = binName(1:2);
electrodes_file = '/mnt/DATA/Clara/ymaze/selected_electrodes2.csv';
electrodes = readtable(electrodes_file);
channelList = electrodes(strcmp(electrodes.animal, animal_code),:).channel;
if isempty(channelList)
    channelList=1:nChans;
end

fs = 600;
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
%dataArray = filter50Hz(dataArray, fs);

freqrange = 1:2:50;
for i = 1:numel(channelList)
    channel = channelList(i);
    figure('Name', [binName '-channel-' num2str(channel)]);
    % Raw signal
    subplot(3,1,1);
    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints, dataArray(channel,:));
    xstd = std(dataArray(channel,:));
    draw_keypoints(time_mouse_arrived, [-5*xstd, 5*xstd], lengthSeconds, secondOffset);
    
    % Spectrogram signal
    [wt, wfreqs]=cwt(dataArray(channel,:), 'amor', fs);
    wt_pow = abs(wt).^2;
    low_freqs = find(wfreqs <= 45);
    high_freqs = find(wfreqs >= 60 & wfreqs < 200);    
    
    subplot(3,1,2);
    draw_cwt(wt_pow(high_freqs,:), timepoints, wfreqs(high_freqs));
    draw_keypoints(time_mouse_arrived, [min(wfreqs(high_freqs)), max(wfreqs(high_freqs))], lengthSeconds, secondOffset);
    
    subplot(3,1,3);
    draw_cwt(wt_pow(low_freqs,:), timepoints, wfreqs(low_freqs));
    draw_keypoints(time_mouse_arrived, [min(wfreqs(low_freqs)), max(wfreqs(low_freqs))], lengthSeconds, secondOffset);
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
end

function [] = draw_cwt(cfs,time,freq)
    A = z_score(cfs);
    args = {time, freq, A};
    surf(args{:},'edgecolor','none');
    view(0,90);
    axis xy;
    axis tight;
    
    shading interp; 
    %colormap(parula(256));
    colormap(jet);
    h = colorbar;
    h.Label.String = 'z-score';
    ylabel('Frequency (Hz)');
end

function [] = draw_laser_points(times, lengthSeconds, secondOffset)
    for i = 1:numel(times)
        x = times(i);
        if x < lengthSeconds && x > 0
            line([x, x], [0 20], 'Color', 'r');
        end
    end
end

function [] = draw_keypoints(time_mouse_arrived, ylim, lengthSeconds, secondOffset)
    if isempty(time_mouse_arrived)
        return
    end
    for i = 1:numel(time_mouse_arrived.sec)
        x = time_mouse_arrived.sec(i) - secondOffset;
        if x < lengthSeconds && x > 0
            line([x, x], ylim, 'Color', 'r');
        end
    end
end
