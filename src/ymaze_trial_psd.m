function [ output ] = ymaze_trial_psd(datarootdir, binfile, showspectrograms )
output = struct(...
    'pow_slow', [],...
    'pow_theta', [],...
    'pow_above_theta', [],...
    'pow_slow_gamma', [],...
    'channel_std', []);

if nargin < 3
    showspectrograms = false;
end

binName = binfile.name;
fprintf('Processing file: %s\n', binName);
meta = ReadMeta(binName, binfile.folder);

secondOffset = 1;
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);
fs = 1250;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
nChans = size(dataArray, 1);

dataArray = filter50Hz(dataArray, fs);

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
    'GoalZone', 'AfterReachedReward_10sec', 'AfterMaze', 'BeforeGoal', 'DuringTask', 'Total'};

key_positions_sample = round(time_mouse_arrived.sec * fs);
reached_reward_sample = time_mouse_arrived.sec(end);
after_reward_sample = reached_reward_sample + 10 * fs;
interval_indecies = [...
    key_positions_sample(:), ...
    [key_positions_sample(2:end); after_reward_sample]];
interval_indecies = [ 
    interval_indecies;
    [(after_reward_sample + round(20 * fs)), size(dataArray,2)];
    [key_positions_sample(1) key_positions_sample(end - 1)];
    [key_positions_sample(1) reached_reward_sample]
    [key_positions_sample(1) size(dataArray,2)]];

slow = [0 4];
theta = [5 10];
slow_gamma = [35 45];
above_theta = [11, (fs / 2 - 1)];

pow_slow = zeros(nChans, numel(key_position_names));
pow_theta = zeros(nChans, numel(key_position_names));
pow_above_theta = zeros(nChans, numel(key_position_names));
pow_slow_gamma = zeros(nChans, numel(key_position_names));
for channel = 1:nChans
    for i = 1:size(interval_indecies, 1)
        left_i = ceil(interval_indecies(i,1));
        right_i = ceil(interval_indecies(i,2));
        if right_i - left_i < fs/4
            if i == 1
                fprintf('Problem with tracking positions in file %s\n', tracking_filepath)
                break
            end
            pow_slow(channel, i) = pow_slow(channel, i - 1);
            pow_theta(channel, i) = pow_theta(channel, i - 1);
            pow_above_theta(channel, i) = pow_above_theta(channel, i - 1);
            pow_slow_gamma(channel, i) = pow_slow_gamma(channel, i - 1);
        else
            [pxx, f1] = pwelch(dataArray(channel,left_i:right_i), ...
                floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
            pow_slow(channel, i) = TotalBandPower(f1, pxx, slow);
            pow_theta(channel, i) = TotalBandPower(f1, pxx, theta);
            pow_above_theta(channel, i) = TotalBandPower(f1, pxx, above_theta);
            pow_slow_gamma(channel, i) = TotalBandPower(f1, pxx, slow_gamma);
        end
    end
end

row_names = cellfun(@(x) num2str(x), num2cell(0:31),'UniformOutput',false);

pow_slow = array2table(pow_slow);
pow_slow.Properties.VariableNames = key_position_names;
pow_slow.Properties.RowNames = row_names;

pow_theta = array2table(pow_theta);
pow_theta.Properties.VariableNames = key_position_names;
pow_theta.Properties.RowNames = row_names;

pow_above_theta = array2table(pow_above_theta);
pow_above_theta.Properties.VariableNames = key_position_names;
pow_above_theta.Properties.RowNames = row_names;

pow_slow_gamma = array2table(pow_slow_gamma);
pow_slow_gamma.Properties.VariableNames = key_position_names;
pow_slow_gamma.Properties.RowNames = row_names;

output.pow_slow = pow_slow;
output.pow_theta = pow_theta;
output.pow_above_theta = pow_above_theta;
output.pow_slow_gamma = pow_slow_gamma;

%% Show spectrogram
if showspectrograms
    freqrange = 1:2:50;
    for channel = 1:nChans
        figure;
        spectrogram(dataArray(channel,:), round(fs / 4), round(fs / 8), freqrange, fs, 'yaxis');
        for i = 1:numel(time_mouse_arrived.sec)
            x = time_mouse_arrived.sec(i);
            line([x, x], [0, max(freqrange)], 'Color', 'r');
        end
    end
end

end

function [ pow ] = TotalBandPower(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,band,'psd');
end
