datarootdir = '/media/prez/DATA/Prez/y-maze/2018-02-12';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
channel = 15;

fprintf('Processing file: %s, channel:%d', binName, channel);

tracking_filepath = get_trackingfilepath(datarootdir, binName);
tracking_dat = readtable(tracking_filepath, 'ReadVariableNames', true);
time_mouse_arrived = readTrackingCsv(tracking_filepath, 0);

%% Set start and end for trace
movie_fs = 15;
first_animal_index = find(tracking_dat.pos_x > -1, 1 );
first_animal_frame = tracking_dat.frame(first_animal_index);
last_animal_index = find(tracking_dat.pos_x > -1, 1, 'last');
last_animal_frame = tracking_dat.frame(last_animal_index);
time_start_sec = first_animal_frame / movie_fs;
time_end_sec = last_animal_frame / movie_fs;

%% Read in data
meta = ReadMeta(binName, path);
lengthSeconds = min(time_end_sec + 8, str2double(meta.fileTimeSecs)) - time_start_sec;
dataArray = ReadSGLXData(meta, time_start_sec, lengthSeconds);
fs=800;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
dataArray = filter50Hz(dataArray, fs);

%% Plot trace + SWRs
times = (1:size(dataArray,2)) / fs;
figure;
subplot(4,1,1);
plot(times, dataArray(channel, :));

% Filter for SWRs
passband = [100 250];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);
filtered = filtfilt(b, a, dataArray(channel,:));
[ripples,sd, normalizedSquaredSignal] = MyFindRipples(times', filtered', ...
                             'frequency', fs, ...
                             'thresholds', [2 4 0.02],...
                             'durations', [30 20 300]);
subplot(4,1,2);
if isempty(ripples)
    plotSWR(times, filtered, fs, [], [])
else
    plotSWR(times, filtered, fs, ripples(:,1), ripples(:,3))
end

subplot(4,1,[3 4]);
%% Calculate spectrogram
times = (1:size(dataArray,2)) / fs;
[wt, wfreqs]=cwt(dataArray(channel,:), 'morse', fs);
draw_cwt(wt, times, wfreqs);

%% Show animal positions
hold on;
times = (first_animal_frame:last_animal_frame) / movie_fs - time_start_sec;
plot(times, 200 - tracking_dat.total_percent(first_animal_index:last_animal_index), 'r')
xlabel('Time (sec)');

function draw_cwt(cfs,time,freq)
    args = {time,freq,abs(cfs).^2};
    surf(args{:},'edgecolor','none');
    view(0,90);
    axis tight;
    shading interp; colormap(parula(256));
    h = colorbar;
    h.Label.String = 'Power';
    xlabel('Time'); ylabel('Hz');
    ylim([0 250]);
end