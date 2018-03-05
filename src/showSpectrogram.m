datarootdir = '/media/prez/DATA/Prez/N&A_rest';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

meta = ReadMeta(binName, path);

secondOffset = 70;
%lengthSeconds = min(50, str2double(meta.fileTimeSecs) - secondOffset);
lengthSeconds = 20;

nChans=32;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 1250;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

dataArray = filter50Hz(dataArray, fs);

tracking_filepath = get_trackingfilepath(datarootdir, binName);
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);

freqrange = 1:2:50;
for channel = 1:nChans
    figure;
    % Raw signal
    subplot(2,1,1);
    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints, dataArray(channel,:));
    xstd = std(dataArray(channel,:));
    draw_keypoints(time_mouse_arrived, [-5*xstd, 5*xstd], lengthSeconds);
    
    % Spectrogram signal
    subplot(2,1,2);
    spectrogram(dataArray(channel,:), round(fs / 4), round(fs / 8), freqrange, fs, 'yaxis');
    draw_keypoints(time_mouse_arrived, [0, max(freqrange)], lengthSeconds);
    b = gca; legend(b,'off');
end

function [] = draw_keypoints(time_mouse_arrived, ylim, lengthSeconds)
    if isempty(time_mouse_arrived)
        return
    end
    for i = 1:numel(time_mouse_arrived.sec)
        x = time_mouse_arrived.sec(i);
        if x < lengthSeconds && x > 0
            line([x, x], ylim, 'Color', 'r');
        end
    end
end