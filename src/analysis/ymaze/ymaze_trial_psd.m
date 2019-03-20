function [ output ] = ymaze_trial_psd(datarootdir, binfile, channels, showspectrograms )
output = struct(...
    'pow',table(),...
    'channel_std', [],...
    'all_psd_xx', [],...
    'wide_psd_xx',[],...
    'psd_all_bands',[],...
    'channel',[]);

if nargin < 4
    showspectrograms = false;
end

binName = binfile.name;
fprintf('Processing file: %s\n', binName);
meta = ReadMeta(binName, binfile.folder);

secondOffset = 0;
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);
fs = 525;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

%dataArray = filter50Hz(dataArray, fs);
nChans = numel(channels);

%% Find noisy channels
output.channel_std = std(dataArray, [], 2);

%% Read in trial info
tracking_filepath = get_trackingfilepath(datarootdir, binName);
time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);
if isempty(time_mouse_arrived)
    return
end

if any(channels > size(dataArray,1))
    return
end

%% Run Power analysis
% Power Spectral analysis

key_position_names = { 'StartZone',  'FirstArm', 'Junction', 'SecondArm', ...
    'GoalZone', 'AfterReachedReward_10sec', 'AfterConsumedReward_10_sec', ...
    'BeforeGoalZone', 'DuringStim', 'DuringMaze', 'HomeCageLast10sec', 'Total'};

key_positions_sample = round(time_mouse_arrived.sec * fs);
key_positions_sample(1) = max(key_positions_sample(1), 1);
if key_positions_sample(1) > size(dataArray, 2)
    return
end

consuming_reward_sec = 10;
after_reward_sample = time_mouse_arrived.sec(end) + consuming_reward_sec * fs;
interval_indecies = [...
    key_positions_sample(:), ... 
    [key_positions_sample(2:end); after_reward_sample]];  % AfterReachedReward_10sec
interval_distances = [time_mouse_arrived.percent(2:end) / 2; 100] - ...
    (time_mouse_arrived.percent(:) / 2);
recording_end_sample = size(dataArray,2);
interval_indecies = [
    interval_indecies;
    [(after_reward_sample + round(10 * fs)), (after_reward_sample + round(20 * fs))];  % AfterConsumerdReward_10_sec
    [key_positions_sample(1) key_positions_sample(end - 1)];  % BeforeGoalZone
    [key_positions_sample(end-1) recording_end_sample]; % DuringStim
    [key_positions_sample(1) key_positions_sample(end)];  % DuringMaze
    [max(recording_end_sample - 10 * fs, key_positions_sample(end-1)) recording_end_sample];  % StimLast10sec
    [key_positions_sample(1) recording_end_sample]  % Total
];
interval_distances = [
    interval_distances;
    0;
    80;
    20;
    100;
    100;
    0;
];


slow = [3 8];
theta = [8 12.5];
slow_gamma = [25 35];
med_gamma = [60 80];
fast_gamma = [80 150];
above_theta = [13, 20];

pow_slow = zeros(nChans, numel(key_position_names));
pow_theta = zeros(nChans, numel(key_position_names));
pow_above_theta = zeros(nChans, numel(key_position_names));
pow_slow_gamma = zeros(nChans, numel(key_position_names));
pow_med_gamma = zeros(nChans, numel(key_position_names));
pow_fast_gamma = zeros(nChans, numel(key_position_names));
dom_freq = zeros(nChans, numel(key_position_names));
theta_dom_freq = zeros(nChans, numel(key_position_names));
ripples_freq = zeros(nChans, numel(key_position_names));
velocity = zeros(nChans, numel(key_position_names));
durations = zeros(nChans, numel(key_position_names));

exp_bands = exp(0.7:0.05:5.3);
all_bands = [2:1:15 exp_bands(43:end)];
all_psd_xx = zeros((numel(all_bands) - 1), numel(key_position_names));
all_ripples = [];


for ch_i = 1:nChans
    
    % Filter for SWRs
    psd_xx = dataArray(channels(ch_i),key_positions_sample(1):end);
    passband = [80 250];
    nyquist = fs / 2;
    filterOrder = 4;
    filterRipple = 20;
    [b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);
    filtered = filtfilt(b, a, psd_xx);
    times = (1:size(filtered,2)) / fs;
    [ripples, ~, ~] = MyFindRipples(times', filtered', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 4 0.01],...
                                 'durations', [10 40 350]);
    nripples = size(ripples, 1);
    if nripples > 0
        translated_ripples = ripples;
%         translated_ripples(:,1) = translated_ripples(:,1) - time_mouse_arrived.sec(end - 1);
%         translated_ripples(:,2) = translated_ripples(:,2) - time_mouse_arrived.sec(end - 1);
%         translated_ripples(:,3) = translated_ripples(:,3) - time_mouse_arrived.sec(end - 1);
        rec_dur_sec = times(end);
        translated_ripples(:,1) = translated_ripples(:,1) - rec_dur_sec;
        translated_ripples(:,2) = translated_ripples(:,2) - rec_dur_sec;
        translated_ripples(:,3) = translated_ripples(:,3) - rec_dur_sec;
        ripples_table = array2table(translated_ripples);
        ripples_table.Properties.VariableNames = {'start_sec', 'peak_t', 'end_time', 'peakpow', 'peak_freq'};

        all_ripples = [all_ripples; ripples_table];
    end
    
    for i = 1:size(interval_indecies, 1)
        left_i = ceil(interval_indecies(i,1));
        right_i = ceil(interval_indecies(i,2));
        left_i = max(1, left_i);
        right_i = min(size(dataArray,2), right_i);
        if right_i - left_i < fs/4
            if i == 1
                fprintf('Problem with tracking positions in file %s\n', tracking_filepath)
                break
            end
%             pow_slow(ch_i, i) = pow_slow(ch_i, i - 1);
%             pow_theta(ch_i, i) = pow_theta(ch_i, i - 1);
%             pow_above_theta(ch_i, i) = pow_above_theta(ch_i, i - 1);
%             pow_slow_gamma(ch_i, i) = pow_slow_gamma(ch_i, i - 1);
%             pow_med_gamma(ch_i, i) = pow_med_gamma(ch_i, i - 1);
%             pow_fast_gamma(ch_i, i) = pow_fast_gamma(ch_i, i - 1);
%             dom_freq(ch_i, i) = dom_freq(ch_i, i-1);
%             theta_dom_freq(ch_i, i) = dom_freq(ch_i, i-1);
%             velocity(ch_i, i) = velocity(ch_i, i-1);
            pow_slow(ch_i, i) = NaN;
            pow_theta(ch_i, i) = NaN;
            pow_above_theta(ch_i, i) = NaN;
            pow_slow_gamma(ch_i, i) = NaN;
            pow_med_gamma(ch_i, i) = NaN;
            pow_fast_gamma(ch_i, i) = NaN;
            dom_freq(ch_i, i) = NaN;
            theta_dom_freq(ch_i, i) = NaN;
            velocity(ch_i, i) = NaN;
            durations(ch_i, i) = 0;
            for j=1:(numel(all_bands)-1)
                all_psd_xx(j, i) = all_psd_xx(j, i-1);
            end
        else
            
            if isempty(ripples)
                section_ripples = [];
            else
                section_ripples = ripples(ripples(:,2) >= left_i/fs & ripples(:,2) <= right_i/fs,:);
                section_ripples(:,1) = section_ripples(:,1) + time_mouse_arrived.sec(1);
                section_ripples(:,2) = section_ripples(:,2) + time_mouse_arrived.sec(1);
                section_ripples(:,3) = section_ripples(:,3) + time_mouse_arrived.sec(1);
            end
            section_sec = (right_i - left_i) / fs;
            ripples_freq(ch_i,i) = size(section_ripples,1) / section_sec;

            psd_xx = dataArray(channels(ch_i),left_i:right_i);
            %[pxx, freqs] = pwelch(psd_xx, ...
            %    floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
            [cws, freqs] = cwt(psd_xx, 'morse', fs);
            pxx = median(abs(cws .^ 2), 2);
            freqs = fliplr(freqs')';
            pxx = fliplr(pxx')';

            [~, maxValIndex] = max(pxx);
            dom_freq(ch_i, i) = freqs(maxValIndex);
            theta_dom_freq(ch_i, i) = DomFreqInBand(freqs, pxx, theta);
            velocity(ch_i, i) = interval_distances(i) / section_sec;
            durations(ch_i, i) = section_sec;
            pow_slow(ch_i, i) = TotalBandPower(freqs, pxx, slow);
            pow_theta(ch_i, i) = TotalBandPower(freqs, pxx, theta);
            pow_above_theta(ch_i, i) = TotalBandPower(freqs, pxx, above_theta);
            pow_slow_gamma(ch_i, i) = TotalBandPower(freqs, pxx, slow_gamma);
            pow_med_gamma(ch_i, i) = TotalBandPower(freqs, pxx, med_gamma);
            pow_fast_gamma(ch_i, i) = TotalBandPower(freqs, pxx, fast_gamma);
            
            for j=1:(numel(all_bands)-1)
                all_psd_xx(j, i) = TotalBandPower(freqs, pxx, [all_bands(j) all_bands(j+1)]);
            end
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
key_position_names = { 'StartZone',  'FirstArm', 'Junction', 'SecondArm', ...
    'GoalZone', 'AfterReachedReward_10sec', 'AfterConsumedReward_10_sec', ...
    'BeforeGoalZone', 'DuringStim', 'DuringMaze', 'HomeCageLast10sec', 'Total'};
pow_med_gamma = array2table(pow_med_gamma);
pow_med_gamma.Properties.VariableNames = strcat('pow_med_gamma_', key_position_names);

pow_fast_gamma = array2table(pow_fast_gamma);
pow_fast_gamma.Properties.VariableNames = strcat('pow_fast_gamma_', key_position_names);

dom_freq = array2table(dom_freq);
dom_freq.Properties.VariableNames = strcat('pow_dom_freq_', key_position_names);

theta_dom_freq = array2table(theta_dom_freq);
theta_dom_freq.Properties.VariableNames = strcat('pow_theta_dom_freq_', key_position_names);

ripples_freq = array2table(ripples_freq);
ripples_freq.Properties.VariableNames = strcat('pow_ripples_freq_', key_position_names);

velocity = array2table(velocity);
velocity.Properties.VariableNames = strcat('velocity_', key_position_names);

durations = array2table(durations);
durations.Properties.VariableNames = strcat('duration_', key_position_names);

%% Replace 0's at Start Zone with power from subsequent stage
for j = (numel(key_position_names) - 1) : -1 : 1
    for i = 1:size(all_psd_xx,1)
        if all_psd_xx(i,j) == 0
            all_psd_xx(i,j) = all_psd_xx(i,j+1);
        end
    end
end

%%
output.pow = [pow_slow pow_theta pow_above_theta pow_slow_gamma pow_med_gamma pow_fast_gamma dom_freq theta_dom_freq ripples_freq velocity durations];
output.all_psd_xx = all_psd_xx ./ repmat(all_psd_xx(:,1), 1, numel(key_position_names));
output.psd_all_bands = all_bands;
output.channel = channel_names;
output.ripples = all_ripples;
%output.ntrials = ntrials;

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

function [ dom_freq ] = DomFreqInBand(freqs, pxx, band)
    
    first_freq_index = find(freqs >= band(1), 1, 'first');
    last_freq_index = find(freqs < band(2), 1, 'last');
    if isempty(last_freq_index) || (last_freq_index - first_freq_index < 1) 
        dom_freq = 0;
        return;
    end
    
    [~, maxValIndex] = max(pxx(first_freq_index:last_freq_index));
    dom_freq = freqs(maxValIndex + first_freq_index - 1);
end

function [ pow ] = TotalBandPower(f1, pxx, band)
    filt_band = [max(min(f1), band(1)) min(max(f1), band(2))];
    if filt_band(2) - filt_band(1) < 0.01
        pow = 0;
        return
    end

    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,filt_band,'psd');
end

