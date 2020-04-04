% TODO:
% - classify as quiet, asleep based on EMG
% - calculate based on diode subtracted channels and compare
% - how to choose electrodes?

%% setup results
ripple_std_thr = 6;
use_diode = 1;
selected_channels_only = 1;
is_urethane = 0;
is_after_ymaze = 0;
datarootdir = '/mnt/DATA/chat_ripples/y-maze';
is_ymaze_trial = 1;
secondOffset = 3;
%datarootdir = '/mnt/DATA/chat_ripples/baseline';
%is_ymaze_trial = 0;
%secondOffset = 0;

trials_fpath = [datarootdir filesep 'trials.csv'];
expstable = readtable(trials_fpath, 'ReadVariableNames', true);
expstable.dirname = strtrim(expstable.dirname);
reverse_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels_reversed.csv';
ord_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels.csv';

nexp = size(expstable, 1);
result_table = table();

all_ripples = table();

% Predefined frequencies
all_bands = exp(0.1:0.07:5.35);
min_section_dur_sec = 0.5;

%% process single experiments
entry_i = 1;
for i = 1:nexp
    animal_code = expstable.animal{i};
    state = expstable.state{i};
        
    dateddir = datestr(expstable.date(i), 'yyyy-mm-dd');
    signalpath = [ datarootdir filesep dateddir filesep strtrim(expstable.dirname{i})];
    binfile = dir([ signalpath filesep animal_code '*.bin']);
    if size(binfile, 1) == 0
        warning('No data files found at %s', signalpath);
        continue
    end
    fprintf('Processing file for date=%s file=%s\n', dateddir, binfile.name);
    meta = ReadMeta(binfile.name, binfile.folder);
    channels_file = reverse_channels_file;
    if ismember('reverse_channel_map', expstable.Properties.VariableNames) &&...
            (expstable.reverse_channel_map(i) == 0)
        channels_file = ord_channels_file;
    end
    channelTable = readChannelTable(...
        channels_file, animal_code, meta, selected_channels_only, use_diode);
    
    dataArray = ReadSGLXData(meta, secondOffset, str2double(meta.fileTimeSecs) - secondOffset);
    dataArray = dataArray(channelTable.rec_order,:);

    fs = 1250;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    time = (1:size(dataArray,2)) / fs;
    if use_diode
        [dataArray, channelTable] = subtractDiodeSignal(dataArray, channelTable);
    end
    
    nchans = size(dataArray, 1);

    %dataArray = filter50Hz(dataArray, fs);
    downfs = 625;
    downsampledDataArray = downsample(dataArray', round(fs / downfs))';

    laserChannelIdx = find(strcmp(channelTable.location, 'Laser'));
    emgIdx = find(strcmp(channelTable.location, 'EMG'));
    if is_ymaze_trial == 0
        trialPeriods = extractTrialPeriodsFromLaser(downsampledDataArray,...
            laserChannelIdx, downfs);
    else
        tracking_filepath = [ datarootdir filesep get_trackingfilepath(dateddir, binfile.name)];
        time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);
        if ~isempty(time_mouse_arrived)
            trialPeriods = createYmazeTrialPeriods(downsampledDataArray,...
                                                   time_mouse_arrived,...
                                                   downfs);
            trialPeriods.laserOn = repmat(expstable.laserOn(i), size(trialPeriods,1), 1);
        else 
            trialPeriods = [];
        end
    end
    
    if is_after_ymaze
        trialPeriods = table();
        trialPeriods.starts = 0;
        trialPeriods.ends = int32(str2double(meta.fileTimeSecs) * downfs);
        trialPeriods.laserOn = expstable.laserOn(i);
        trialPeriods.stage_desc = {'after'};
    end
    
    if isempty(trialPeriods)
        warning('No trial periods found in the recording')
        continue
    end
    
    % Filter LFP for SWR detection
    filtered = applyRippleFilter(dataArray, channelTable, fs);

    %% calculate PSD

    % Freq bands during Ymaze trial and Mobility
    slow = [3 7];
    theta = [7 12];
    above_theta = [12 25];
    slow_gamma = [25 45];
    med_gamma = [60 80];
    fast_gamma = [80 130];
    ripple_band = [130 250];
    
    % Urethane specific bounds
    if is_urethane
        slow = [0.1 2];
        theta = [2 5];
        above_theta = [7 20];
        slow_gamma = [20 40];
    else
        if ~is_ymaze_trial % Immobility
            slow = [2 3];
            theta = [5 8];
            above_theta = [10 20];
            slow_gamma = [25 45];
        end
    end

    keep_sample = excludeEMGNoisePeriods(dataArray(emgIdx,:), fs * 0.5);
    keep_sample_fewer = excludeEMGNoisePeriods(dataArray(emgIdx,:), fs * 1);
    for channel = 1:nchans
        if strcmp(channelTable.location{channel}, 'Laser') || ...
                strcmp(channelTable.location{channel}, 'EMG')
            continue
        end
        ripple_detection_signal = GetRippleSignal(filtered(channel, :)', fs);
%         std_estimate_index = find(...
%             strcmp(std_estimates.animal, animal_code) & ...
%             std_estimates.date == expstable.date(i) & ...
%             strcmp(strrep(std_estimates.channel_name, ' ', ''), ...
%                    strrep(channelTable.channel_name{channel}, ' ', '')));
        
        % remove epochs with jumps in the signal
        channel_keep_sample = excludeEMGNoisePeriods(dataArray(channel,:), fs * 0.5);
        keep_sample = intersect(channel_keep_sample, keep_sample);
        %std_estimate = std_estimates.std_estimate(std_estimate_index);
        %std_estimate = median(ripple_detection_signal(keep_sample_fewer) / 0.6745);
        if numel(keep_sample_fewer) < 2
            keep_sample_fewer = keep_sample;
        end
        std_estimate = std(ripple_detection_signal(keep_sample_fewer));
        fprintf('Using std estimate %.8f\n', std_estimate);
        [ripples, sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(channel,:)', ...
                     'frequency', fs, ...
                     'thresholds', [2 ripple_std_thr 0.01],...
                     'durations', [10 20 300],...
                     'std', std_estimate);
        if ~isempty(ripples)
            keep_ripples = ismember(int32(ripples(:,2) * fs), keep_sample);
            if ~all(keep_ripples)
                warning('Rejected %d ripple(s) because of the noise detected in EMG', ...
                    sum(keep_ripples==0));
            end
            ripples = ripples(keep_ripples,:);
        end
        
        %[pxx, freqs] = pwelch(downsampledDataArray, ...
        %    floor(downfs / 4), floor(downfs / 8), floor(downfs / 2), downfs);
        [cws, freqs] = cwt(downsampledDataArray(channel,:), 'amor', downfs);        
        signal_pxx = abs(cws .^ 2);
        freqs = fliplr(freqs')';
        signal_pxx = fliplr(signal_pxx')';
        
        band_power = zeros(numel(all_bands), size(signal_pxx, 2));
        for j=1:(numel(all_bands)-1)
            if all_bands(j) >= freqs(1)
                % Assign power of the first higher frequency
                k = find(freqs >= all_bands(j), 1, 'first');
                band_power(j,:) = signal_pxx(k, :);
            end
        end
        band_power_zscored = zscore(band_power);
        
        for trial_period_i = 1:size(trialPeriods, 1)
            period_start = max(1, int32(trialPeriods.starts(trial_period_i)));
            period_end = min( size(downsampledDataArray, 2), ...
                              int32(trialPeriods.ends(trial_period_i)));
            period_signal = downsampledDataArray(channel, period_start : period_end);
            period_pxx = signal_pxx(:, period_start : period_end);
            sec_start = max(0, double(period_start) / downfs);
            sec_end = double(period_end) / downfs;
            %sec_end = min(sec_end,  size(downsampledDataArray, 2) / downfs);
            sec_length = sec_end - sec_start;

            if sec_length <= min_section_dur_sec
                continue
            end

            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            result_table.channel(entry_i) = channelTable.channel(channel);
            result_table.channelLocation(entry_i) = channelTable.location(channel);
            result_table.channelName(entry_i) = channelTable.channel_name(channel);
            result_table.date(entry_i) = {dateddir};
            result_table.animal(entry_i) = {animal_code};
            result_table.file_name(entry_i) = expstable.dirname(i);
            trial = [expstable.dirname{i} '_' num2str(idivide(int16(trial_period_i-1),3))];
            result_table.trial(entry_i) = {trial};
            result_table.stage_desc(entry_i) = trialPeriods.stage_desc(trial_period_i);
            result_table.state(entry_i) = {state};
            result_table.laserOn(entry_i) = trialPeriods.laserOn(trial_period_i);
            

            slow_adjusted = [ max(slow(1), freqs(1)) slow(2) ];
            result_table.dom_freq(entry_i) = PeakFreq(freqs, period_pxx, [slow_adjusted(1), 45]);

            %result_table.pow_slow(entry_i) = MeanBandPower(freqs, period_pxx, slow_adjusted);
            result_table.pow_slow(entry_i) = CalcBandPower(downfs, period_signal, slow_adjusted);
            %result_table.pow_theta(entry_i) = MeanBandPower(freqs, period_pxx, theta);
            result_table.pow_theta(entry_i) = CalcBandPower(downfs, period_signal, theta);
            result_table.peak_theta(entry_i) = PeakFreq(freqs, period_pxx, theta);
            %result_table.pow_slow_gamma(entry_i) = MeanBandPower(freqs, period_pxx, slow_gamma);
            result_table.pow_slow_gamma(entry_i) = CalcBandPower(downfs, period_signal, slow_gamma);
            result_table.peak_slow_gamma(entry_i) = PeakFreq(freqs, period_pxx, slow_gamma);
            %result_table.pow_med_gamma(entry_i) = MeanBandPower(freqs, period_pxx, med_gamma);
            result_table.pow_med_gamma(entry_i) = CalcBandPower(downfs, period_signal, med_gamma);
            %result_table.pow_fast_gamma(entry_i) = MeanBandPower(freqs, period_pxx, fast_gamma);
            result_table.pow_fast_gamma(entry_i) = CalcBandPower(downfs, period_signal, fast_gamma);
            %result_table.pow_above_theta(entry_i) = MeanBandPower(freqs, period_pxx, above_theta);
            result_table.pow_above_theta(entry_i) = CalcBandPower(downfs, period_signal, above_theta);
            %result_table.pow_ripple_band(entry_i) = MeanBandPower(freqs, period_pxx, ripple_band);
            result_table.pow_ripple_band(entry_i) = CalcBandPower(downfs, period_signal, ripple_band);
            for j=1:(numel(all_bands)-1)
                result_table.all_psd_xx(entry_i, j) = ...
                    mean(band_power(j, period_start : period_end));
                result_table.all_psd_z(entry_i, j) = ...
                    mean(band_power_zscored(j, period_start : period_end));
            end

             section_ripples = [];
             if ~isempty(ripples)
                section_ripples = ripples(ripples(:,2) >= sec_start & ripples(:,2) <= sec_end,:);
             end
             nripples = size(section_ripples, 1);
             if ~isempty(section_ripples)
                result_table.has_ripples(entry_i) = 1;
                result_table.nripples(entry_i) = size(section_ripples,1);
                result_table.swr_incidence(entry_i) = size(section_ripples,1) / sec_length;
                result_table.ripple_peakpow(entry_i) = mean(section_ripples(:,4));
                result_table.ripple_dur(entry_i) = mean(section_ripples(:,3) - section_ripples(:,1));
                result_table.ripple_freq(entry_i) = mean(section_ripples(:,5));
               
                % Align ripple timings to the start of the trial
                section_ripples(:, 1:3) = section_ripples(:, 1:3) - sec_start;
                ripples_table = array2table(section_ripples);
                ripples_table.Properties.VariableNames = {'start_sec', 'peak_t', 'end_time', 'peakpow', 'peak_freq'};
                ripples_table.abs_peak_t = ripples_table.peak_t + sec_start;
                ripples_table.date = repmat(dateddir, nripples, 1);
                ripples_table.animal = repmat(animal_code, nripples, 1);
                ripples_table.file_name = repmat(expstable.dirname(i), nripples, 1);
                ripples_table.trial = repmat({trial}, nripples, 1);
                ripples_table.stage_desc = repmat({trialPeriods.stage_desc{trial_period_i}}, nripples, 1);
                ripples_table.stage_dur_sec = repmat(sec_length, nripples, 1);
                ripples_table.channel = repmat(channelTable.channel(channel), nripples, 1);
                ripples_table.channelLocation = repmat(channelTable.location(channel), nripples, 1);
                ripples_table.channelName = repmat(channelTable.channel_name(channel), nripples, 1);
                ripples_table.laserOn = repmat(result_table.laserOn(entry_i), nripples, 1);
                ripples_table.state = repmat({state}, nripples, 1);
 
                all_ripples = [all_ripples; ripples_table];
             end
             if ~isempty(ripples)
                result_table.swr_powerhz(entry_i) = mean(abs(normalizedSquaredSignal));
             end
            
             entry_i = entry_i + 1;
        end
       
    end
end

filename_infix = '';
if use_diode
    filename_infix = '_diode';
end

if is_after_ymaze
    filename_infix = [filename_infix '_after'];
end
if ~selected_channels_only
    filename_infix = [filename_infix '_all'];
end

%filename_infix = [filename_infix '_single_std'];
outfile_suffix = [filename_infix '_th' num2str(ripple_std_thr) '.csv'];
writetable(result_table, [datarootdir filesep 'psd_table' outfile_suffix]);
writetable(all_ripples, [datarootdir filesep 'ripples' outfile_suffix]);
            
function [ pow ] = TotalBandPower(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,band,'psd');
end

function [ pow ] = CalcBandPower(fs, signal, band)
    pow = bandpower(signal, fs, band);
    %band_freqs = find(f1 >= band(1) & f1 <= band(2));
    %pow = trapz(f1(band_freqs), pxx(band_freqs)) / abs(band(2) - band(1));
end

function [ pow ] = TotalBandPower2(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    band_freqs = find(f1 >= band(1) & f1 <= band(2));
    pow = -trapz(f1(band_freqs), log10(pxx(band_freqs))) / (band(2) - band(1));
end

function [ freq ] = PeakFreq(f1, pxx, band)
    pxx = mean(pxx, 2);
    band_freqs_index = find(f1 >= band(1) & f1 <= band(2));
    [~, maxValIndex] = max(pxx(band_freqs_index));
    freqs = f1(band_freqs_index);
    freq = freqs(maxValIndex);
end

