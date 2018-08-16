datarootdir = 'F:\LFP_data\y-maze\2018-02-13';
path = [datarootdir filesep 'signal'];
[binName, path] = uigetfile('*.bin', 'LFP file', path);
channel = 17;

fprintf('Processing file: %s, channel:%d', binName, channel);

tracking_filepath = get_trackingfilepath(datarootdir, binName);
tracking_dat = readtable(tracking_filepath, 'ReadVariableNames', true);
time_mouse_arrived = readTrackingCsv(tracking_filepath, 0);

%% Set start and end for trace
movie_fs = 15;
first_animal_index = find(tracking_dat.pos_x > -1, 1 ) + round(0.6 * movie_fs);
first_animal_frame = tracking_dat.frame(first_animal_index);
last_animal_index = find(tracking_dat.pos_x > -1, 1, 'last');
last_animal_frame = tracking_dat.frame(last_animal_index);
time_start_sec = first_animal_frame / movie_fs + 3;
time_end_sec = last_animal_frame / movie_fs;

%% Read in data
meta = ReadMeta(binName, path);
seconds_after_goal = 10;
lengthSeconds = min(time_end_sec + seconds_after_goal, str2double(meta.fileTimeSecs)) - time_start_sec;
dataArray = ReadSGLXData(meta, time_start_sec, lengthSeconds);
fs=800;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
dataArray = filter50Hz(dataArray, fs);

figure;
%% Show animal positions
subplot(4,1,1);
times = (first_animal_frame:last_animal_frame) / movie_fs - time_start_sec;
plot(times, tracking_dat.total_percent(first_animal_index:last_animal_index) / 2, 'r')
xlim([0, lengthSeconds]);
ylim([0,100]);
tracking_filepath = get_trackingfilepath(datarootdir, binName);
hold on;
draw_keypoints(time_mouse_arrived, [0 200], lengthSeconds, time_start_sec);

%% Plot trace
times = (1:size(dataArray,2)) / fs;
subplot(4,1,2);
plot(times, dataArray(channel, :));


subplot(4,1,[3 4]);
%% Calculate spectrogram
times = (1:size(dataArray,2)) / fs;
[cfs, wfreqs]=cwt(dataArray(channel,:), 'morse', fs);
draw_cwt(cfs, times, wfreqs);
xlabel('Time (sec)');


%% Calculate PSD at goal zone vs start zone
figure('Name', 'Start vs at reward');
start_index = max(1, time_mouse_arrived.sec(end-2) * fs - time_start_sec * fs);
end_index = time_mouse_arrived.sec(end-1) * fs;
pxx = mean(abs(cfs(:,start_index:end_index) .^ 2), 2);
semilogx(wfreqs,log10(pxx));
hold on;

start_index = (time_mouse_arrived.sec(end-1) - time_start_sec) * fs;
end_index = (time_mouse_arrived.sec(end) - time_start_sec) * fs;
pxx = mean(abs(cfs(:,start_index:end_index) .^ 2), 2);
semilogx(wfreqs,log10(pxx));
xlim([1,200]);



%% Functions
function plot_band(times, pxx, freqs, band) 
    first_freq_index = find(freqs <= band(2), 1, 'first');
    last_freq_index = find(freqs >= band(1), 1, 'last');
    band_pow = sum(pxx(first_freq_index:last_freq_index,:));
    
    % smooth the signal pow
    w = gausswin(1000);
    smoothed_band_pow = filter(w, 1, band_pow);
    plot(times, smoothed_band_pow)
end


function [] = draw_cwt(cfs,time,freq)
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