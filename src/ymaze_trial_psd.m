function [ output ] = ymaze_trial_psd(datarootdir, binfile, channels, showspectrograms )
output = struct(...
    'pow',table(),...
    'channel_std', []);

if nargin < 4
    showspectrograms = false;
end

binName = binfile.name;
fprintf('Processing file: %s\n', binName);
meta = ReadMeta(binName, binfile.folder);

secondOffset = 1;
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);
fs = 525;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

dataArray = filter50Hz(dataArray, fs);
nChans = numel(channels);

%% Find noisy channels
output.channel_std = std(dataArray, [], 2);

%% Read in trial info
tracking_filepath = get_trackingfilepath(datarootdir, binName);
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);
if isempty(time_mouse_arrived)
    return
end

%% Run Power analysis
% Power Spectral analysis

key_position_names = { 'StartZone',  'FirstArm', 'Junction', 'SecondArm', ...
    'GoalZone', 'AfterReachedReward_10sec', 'AfterConsumedReward_10_sec', ...
    'BeforeGoal', 'DuringMaze', 'HomeCageLast10sec', 'Total'};

key_positions_sample = round(time_mouse_arrived.sec * fs);
reached_reward_sample = time_mouse_arrived.sec(end);
consuming_reward_sec = 10;
after_reward_sample = time_mouse_arrived.sec(end-1) + consuming_reward_sec * fs;
interval_indecies = [...
    key_positions_sample(:), ...
    [key_positions_sample(2:end); after_reward_sample]];
interval_distances = [time_mouse_arrived.percent(2:end) / 2; 100] - ...
    (time_mouse_arrived.percent(:) / 2);
recording_end_sample = size(dataArray,2);
interval_indecies = [
    interval_indecies;
    [(after_reward_sample + round(10 * fs)), size(dataArray,2)];
    [key_positions_sample(1) key_positions_sample(end - 1)];
    [key_positions_sample(1) reached_reward_sample]
    [key_positions_sample(1) recording_end_sample]
    (recording_end_sample - 10 * fs) recording_end_sample];
interval_distances = [
    interval_distances;
    0;
    80;
    100;
    100;
    0;
];

slow = [0.8 5];
theta = [6 10];
slow_gamma = [30 45];
med_gamma = [60 120];
fast_gamma = [100 200];
above_theta = [10.5, 200];

pow_slow = zeros(nChans, numel(key_position_names));
pow_theta = zeros(nChans, numel(key_position_names));
pow_above_theta = zeros(nChans, numel(key_position_names));
pow_slow_gamma = zeros(nChans, numel(key_position_names));
pow_med_gamma = zeros(nChans, numel(key_position_names));
pow_fast_gamma = zeros(nChans, numel(key_position_names));
dom_freq = zeros(nChans, numel(key_position_names));
ripples_freq = zeros(nChans, numel(key_position_names));
velocity = zeros(nChans, numel(key_position_names));

for ch_i = 1:nChans
    for i = 1:size(interval_indecies, 1)
        left_i = ceil(interval_indecies(i,1));
        right_i = ceil(interval_indecies(i,2));
        if right_i - left_i < fs/4
            if i == 1
                fprintf('Problem with tracking positions in file %s\n', tracking_filepath)
                break
            end
            pow_slow(ch_i, i) = pow_slow(ch_i, i - 1);
            pow_theta(ch_i, i) = pow_theta(ch_i, i - 1);
            pow_above_theta(ch_i, i) = pow_above_theta(ch_i, i - 1);
            pow_slow_gamma(ch_i, i) = pow_slow_gamma(ch_i, i - 1);
            pow_med_gamma(ch_i, i) = pow_med_gamma(ch_i, i - 1);
            pow_fast_gamma(ch_i, i) = pow_fast_gamma(ch_i, i - 1);
            dom_freq(ch_i, i) = dom_freq(ch_i, i-1);
            velocity(ch_i, i) = velocity(ch_i, i-1);
        else
            psd_xx = dataArray(channels(ch_i),left_i:right_i);

            % Filter for SWRs
            passband = [100 250];
            nyquist = fs / 2;
            filterOrder = 4;
            filterRipple = 20;
            [b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);
            filtered = filtfilt(b, a, psd_xx);
            times = (1:size(filtered,2)) / fs;
            [ripples,sd, normalizedSquaredSignal] = MyFindRipples(times', filtered', ...
                                         'frequency', fs, ...
                                         'thresholds', [2 4 0.02],...
                                         'durations', [30 20 300]);
            ripples_freq(ch_i,i) = size(ripples,1) / (times(end) - times(1));

            %[pxx, freqs] = pwelch(psd_xx, ...
            %    floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
            [cws, freqs] = cwt(psd_xx, 'morse', fs);
            pxx = mean(abs(cws .^ 2), 2);
            freqs = fliplr(freqs')';
            pxx = fliplr(pxx')';

            [~, maxValIndex] = max(pxx);
            dom_freq(ch_i, i) = freqs(maxValIndex);

            velocity(ch_i, i) = interval_distances(i) / (times(end) - times(1));
            pow_slow(ch_i, i) = TotalBandPower(freqs, pxx, slow);
            pow_theta(ch_i, i) = TotalBandPower(freqs, pxx, theta);
            pow_above_theta(ch_i, i) = TotalBandPower(freqs, pxx, above_theta);
            pow_slow_gamma(ch_i, i) = TotalBandPower(freqs, pxx, slow_gamma);
            pow_med_gamma(ch_i, i) = TotalBandPower(freqs, pxx, med_gamma);
            pow_fast_gamma(ch_i, i) = TotalBandPower(freqs, pxx, fast_gamma);
        end
    end
end

channel_names = cellfun(@(x) num2str(x), num2cell(channels),'UniformOutput',false);

pow_slow = array2table(pow_slow);
pow_slow.Properties.VariableNames = strcat('pow_slow_', key_position_names);

pow_theta = array2table(pow_theta);
pow_theta.Properties.VariableNames = strcat('pow_theta_', key_position_names);

pow_above_theta = array2table(pow_above_theta);
pow_above_theta.Properties.VariableNames = strcat('pow_above_theta_', key_position_names);

pow_slow_gamma = array2table(pow_slow_gamma);
pow_slow_gamma.Properties.VariableNames = strcat('pow_slow_gamma_', key_position_names);

pow_med_gamma = array2table(pow_med_gamma);
pow_med_gamma.Properties.VariableNames = strcat('pow_med_gamma_', key_position_names);

pow_fast_gamma = array2table(pow_fast_gamma);
pow_fast_gamma.Properties.VariableNames = strcat('pow_fast_gamma_', key_position_names);

dom_freq = array2table(dom_freq);
dom_freq.Properties.VariableNames = strcat('pow_dom_freq_', key_position_names);

ripples_freq = array2table(ripples_freq);
ripples_freq.Properties.VariableNames = strcat('pow_ripples_freq_', key_position_names);

velocity = array2table(velocity);
velocity.Properties.VariableNames = strcat('velocity_', key_position_names);

output.pow = [pow_slow pow_theta pow_above_theta pow_slow_gamma pow_med_gamma pow_fast_gamma dom_freq ripples_freq velocity];
output.channel = channel_names;

%% Show spectrogram
if showspectrograms
    freqrange = 1:2:50;
    for ch_i = 1:nChans
        figure;
        spectrogram(dataArray(ch_i,:), round(fs / 4), round(fs / 8), freqrange, fs, 'yaxis');
        for i = 1:numel(time_mouse_arrived.sec)
            x = time_mouse_arrived.sec(i);
            line([x, x], [0, max(freqrange)], 'Color', 'r');
        end
    end
end

end

function [ pow ] = TotalBandPower(f1, pxx, band)
    filt_band = [max(min(f1), band(1)) min(max(f1), band(2))];
    if filt_band(2) - filt_band(1) < 0.5
        pow = 0;
        return
    end

    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,filt_band,'psd');
end

