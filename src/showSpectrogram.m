datarootdir = '/media/prez/DATA/Prez/N&A_rest';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

meta = ReadMeta(binName, path);

secondOffset = 0;
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
%lengthSeconds = 40;

nChans=32;
channelList = [16 2 1];

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 120;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

dataArray = filter50Hz(dataArray, fs);

tracking_filepath = get_trackingfilepath(datarootdir, binName);
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);

statefilepath = [datarootdir filesep 'state.csv'];
statetable = readtable(statefilepath, 'ReadVariableNames', true);
binNameParts = strsplit(binName, '_g0');
expname = binNameParts{1};
expstates_index = find(...
    cellfun(@(x) strcmp(x,expname), statetable.experiment, 'UniformOutput',1));


times = (1:size(dataArray,2)) / fs + secondOffset;
%freqrange = 1:2:50;
for i = 1:size(channelList,2)
    channel = channelList(i);
    figure('Name', [binName '-channel-' num2str(channel)]);
    % Raw signal
    subplot(2,1,1);
    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints, dataArray(channel,:));
    xstd = std(dataArray(channel,:));
    draw_keypoints(time_mouse_arrived, [-5*xstd, 5*xstd], lengthSeconds);
    
    % Spectrogram signal
    subplot(2,1,2);
    [wt, wfreqs]=cwt(dataArray(channel,:), 'morse', fs);
    draw_cwt(wt, times, wfreqs);
    %spectrogram(dataArray(channel,:), round(fs / 4), round(fs / 8), freqrange, fs, 'yaxis');
    %draw_keypoints(time_mouse_arrived, [0, max(freqrange)], lengthSeconds);
    draw_laser_points(statetable.time_sec(expstates_index), lengthSeconds);
    b = gca; legend(b,'off');
end

function draw_cwt(cfs,time,freq)
    args = {time,freq,abs(cfs).^2};
    surf(args{:},'edgecolor','none');
    view(0,90);
    axis tight;
    shading interp; colormap(parula(256));
    h = colorbar;
    h.Label.String = 'Power';
    xlabel('Time'); ylabel('Hz');
    ylim([0 40]);
end

function [] = draw_laser_points(times, lengthSeconds)
    for i = 1:numel(times)
        x = times(i);
        if x < lengthSeconds && x > 0
            line([x, x], [0 20], 'Color', 'r');
        end
    end
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