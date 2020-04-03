is_ymaze_trial = 1;
animal_code = 'OS';

if is_ymaze_trial
    datarootdir = '/mnt/DATA/chat_ripples/y-maze';
    secondOffset = 2;
    trial_date = datetime('2019-11-13', 'InputFormat', 'yyyy-MM-dd');
    laserOn = 1;
else
    datarootdir = '/mnt/DATA/chat_ripples/baseline';
    secondOffset = 0;
    laserOn = 0;
end

trials_fpath = [datarootdir filesep 'trials.csv'];
expstable = readtable(trials_fpath, 'ReadVariableNames', true);
expstable.dirname = strtrim(expstable.dirname);
reverse_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels_reversed_same_side.csv';
ord_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels.csv';
nexp = size(expstable, 1);

if is_ymaze_trial
    daytable = expstable(...
        expstable.date >= trial_date & ...
        strcmp(expstable.animal, repmat({animal_code}, nexp, 1)) &...
        expstable.laserOn == laserOn & expstable.correct == 1, :);
else
    daytable = expstable(...
         strcmp(expstable.animal, repmat({animal_code}, nexp, 1)), :);
    %daytable = expstable;
end

result_table = table();

all_ripples = table();

all_bands = exp(0.7:0.05:5.55);

if is_ymaze_trial
    trial_period_lengths = [10 20 50 50];
else
    trial_period_lengths = [20 20 20];
end

% Matrix with spectrogram values indexed: channel location x 
% freq band index x common timestamp
band_zscore = zeros(2, numel(all_bands)-1, sum(trial_period_lengths));
% Count of trials: trial period x common timestamp
ntrials = zeros(1, sum(trial_period_lengths));

band_coherence = zeros(numel(all_bands)-1, sum(trial_period_lengths));


channelLocs = {'CA1', 'CA3'};
%% process single experiments
entry_i = 1;
for i = 1:size(daytable, 1)
    animal_code = daytable.animal{i};
    date_str = datestr(daytable.date(i), 'yyyy-mm-dd');
    signalpath = [ datarootdir filesep date_str filesep strtrim(daytable.dirname{i})];
    binfile = dir([ signalpath filesep animal_code '*.bin']);
    if size(binfile, 1) == 0
        warning('No data files found at %s', signalpath);
        continue
    end
    fprintf('Processing file for date=%s file=%s\n', date_str, binfile.name);
    meta = ReadMeta(binfile.name, binfile.folder);
    channels_file = reverse_channels_file;
    if ismember('reverse_channel_map', daytable.Properties.VariableNames) &&...
            (daytable.reverse_channel_map(i) == 0)
        channels_file = ord_channels_file;
    end
    channelTable = readChannelTable(...
        channels_file, animal_code, meta, 1, 0);
    laserChannelIdx = find(strcmp(channelTable.location, 'Laser'));
    emgIdx = find(strcmp(channelTable.location, 'EMG'));
    
    dataArray = ReadSGLXData(meta, secondOffset, str2double(meta.fileTimeSecs) - secondOffset);
    dataArray = dataArray(channelTable.rec_order,:);

    fs = 625;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    time = (1:size(dataArray,2)) / fs;
    nchans = size(dataArray, 1);    
    
    if is_ymaze_trial == 1
        tracking_filepath = [ datarootdir filesep get_trackingfilepath(date_str, binfile.name)];
        time_mouse_arrived = readTrackingCsv(tracking_filepath, secondOffset);
        if ~isempty(time_mouse_arrived)
            trialPeriods = createYmazeTrialPeriods(dataArray,...
                                                   time_mouse_arrived,...
                                                   fs);
            trialPeriods.laserOn = repmat(expstable.laserOn(i), size(trialPeriods,1), 1);
            trialPeriods = trialPeriods(...
                strcmp(trialPeriods.stage_desc, 'StartZone') | ...
                strcmp(trialPeriods.stage_desc, 'MazeStem') | ...
                strcmp(trialPeriods.stage_desc, 'StimLast10sec') | ...
                strcmp(trialPeriods.stage_desc, 'GoalZone20sec'), :);
            trialPeriods = sortrows(trialPeriods, 'starts');
            trialPeriods.trial_ordinal = ones(size(trialPeriods,1), 1);
        else 
            continue;
        end
    else
        trialPeriods = extractTrialPeriodsFromLaser(dataArray, laserChannelIdx, fs);
        %trialPeriods = trialPeriods(~strcmp(trialPeriods.stage_desc, 'after_stim'), :);
        sortrows(trialPeriods, 'starts');
        trialPeriods.trial_ordinal = (repelem(1:(size(trialPeriods,1)/3), 1, 3))';
    end

    % Calculate coherence for the predefined frequencies
    ca_channels = find(startsWith(channelTable.location, 'CA'));
    
    for trial_ordinal = (unique(trialPeriods.trial_ordinal))'
        trialSections = trialPeriods(trialPeriods.trial_ordinal == trial_ordinal, :);
        
        section_start = max(1, min(trialSections.starts));
        section_end = min(max(trialSections.ends), size(dataArray,2));
        section_signal = dataArray(:, section_start : section_end);
        
        [wcoh, wcs, freqs] = wcoherence(section_signal(ca_channels(1),:), ...
            section_signal(ca_channels(2),:), fs);
        freqs = fliplr(freqs')';
        wcoh = fliplr(wcoh')';
        band_coh = zeros(numel(all_bands)-1, size(wcoh, 2));
        for j=1:(numel(all_bands)-1)
            % Assign power of the first higher frequency
            k = find(freqs >= all_bands(j), 1, 'first');
            band_coh(j,:) = wcoh(k, :);
        end
        emgIdx = find(strcmp(channelTable.location, 'EMG'));
        periods_coherence = zeros(size(band_coherence, 1), size(band_coherence, 2));

        %dataArray = dataArray(:, keep_sample);
        for chan_i = 1:size(channelTable, 1)
            loc = channelTable.location(chan_i);

            if strcmp(loc, 'EMG') || strcmp(loc, 'Laser')
                continue
            end
            chan_output_i = find(strcmp(loc, channelLocs));

            [wt, wfreqs] = cwt( section_signal(chan_i,:), 'morse', fs, ...
                'ExtendSignal', true, ...
                'VoicePerOctave', 30, 'WaveletParameters', [3 120]);
            pxx = abs(wt).^2;
            wfreqs = fliplr(wfreqs')';
            pxx = fliplr(pxx')';

            band_power = zeros(numel(all_bands)-1, size(pxx,2));
            for j=1:(numel(all_bands)-1)
                if all_bands(j) >= wfreqs(1)
                    band_power(j, :) = bandpower(pxx, ...
                        wfreqs, [all_bands(j) all_bands(j+1)], 'psd');
                end
            end

            periods_bandpower = zeros(size(band_zscore, 2), size(band_zscore, 3));
            for section_i = 1:size(trialSections, 1)
                period_start_i = max(1, trialSections.starts(section_i) - section_start + 1);
                period_end_i = min(trialSections.ends(section_i) - section_start + 1,...
                               size(section_signal,2));
                if period_end_i - period_start_i <= 0
                    continue
                end

                trial_period_offsets = cumsum([0 trial_period_lengths]);
                offset_i = trial_period_offsets(section_i);
                period_len = trial_period_lengths(section_i);

                keep_sample = excludeEMGNoisePeriods(dataArray(emgIdx,period_start_i:period_end_i), fs * 0.2);
                period_pxx = band_power(:, period_start_i:period_end_i);
                period_indecies = (offset_i + 1):(offset_i + period_len);
                periods_bandpower(:, period_indecies) = ...
                    matScale(period_pxx, period_len);

                period_coh = band_coh(:, period_start_i:period_end_i);
                periods_coherence(:, period_indecies) = ...
                    matScale(period_coh, period_len);

                ntrials(period_indecies) = ntrials(period_indecies) + 1;
            end
            band_zscored_trial = zscore(periods_bandpower);
            %band_zscored_trial = movmean(band_zscored_trial, 3, 2);
            band_zscore(chan_output_i, :, :) = ...
                squeeze(band_zscore(chan_output_i, :, :)) + band_zscored_trial;

        end
        band_coherence = band_coherence + periods_coherence;
    end
end

%% Plot spectrograms
xintercepts = [1 cumsum(trial_period_lengths)];
if is_ymaze_trial
    xlabels = {'Start Zone','Maze Stem','Goal Zone', 'Last 10 sec'};
else
    xlabels = {'laser OFF', 'laser ON', 'laser OFF'};
end
    
file_prefix = 'sleep_';
if is_ymaze_trial
    file_prefix = 'ymaze_';
end
file_prefix = [file_prefix animal_code '_'];

for chan_i = [1 2]
    f = figure('name', ['Channel=' channelLocs{chan_i} ],...
               'Position', [400 700 800 300]);
    zscore_mean = squeeze(band_zscore(chan_i, :, :)) ./ ...
        squeeze((ntrials / 2)); % Each trial counted twice: once per channel
    %zscore_mean = movmean(zscore_mean, 5, 2);
    zscore_mean = apply_ylim(zscore_mean);
    plotSpectrogram(zscore_mean, 1:sum(trial_period_lengths), ...
        all_bands(2:end), xintercepts, xlabels);
    saveas(f, ['/home/prez/tmp/ymaze_chat_x_ai32/' file_prefix 'spect' channelLocs{chan_i} '.svg'])
end

%% Plot coherence
% Each trial counted twice: once per channel
scaled_band_coherence = zscore(band_coherence ./ (ntrials / 2));
f = figure('name', 'Wavelet Coherence',...
           'Position', [400 700 800 300]);

plotSpectrogram(scaled_band_coherence, 1:sum(trial_period_lengths), ...
    all_bands(2:end), xintercepts, xlabels);
saveas(f, ['/home/prez/tmp/ymaze_chat_x_ai32/' file_prefix '_coherence_laser' num2str(laserOn) '.svg'])

%% Utilify functions
function [scaledM] = matScale(M, required_width)
    % Stretch the trial period length
    mwidth = size(M, 2);
    if mwidth < required_width
        M = repelem(M, 1, ceil(required_width / mwidth));
    end

    % Compress the trial period
    blockSize = [1, floor(mwidth / required_width)];
    meanFilterFun = @(blockStruct) mean2(blockStruct.data(:));
    scaledM = blockproc(M, blockSize, meanFilterFun);
    scaledM = scaledM(:, 1:required_width);
end

function [] = plotSpectrogram(wt_pow, time, wfreqs, xintercepts, xlabels)
    low_freqs = find(wfreqs <= 15);
    high_freqs = find(wfreqs < 250 & wfreqs > 15);
    
    % Workaround to force saving as svg with all the elements
    set(0, 'DefaultFigureRenderer', 'painters');
    
    subplot(3,1,1:2);
    draw_cwt(wt_pow(high_freqs,:), time, wfreqs(high_freqs));
    h = colorbar;
    h.Label.String = 'z-score';
    
    xticks(xintercepts);

    subplot(3,1,3);
    draw_cwt(wt_pow(low_freqs,:), time, wfreqs(low_freqs));
    h = colorbar;
    h.Label.String = 'z-score';
    %draw_vlines(xintercepts, [min(wfreqs(low_freqs)) min(wfreqs(low_freqs)) + 0.6]);
    ax = gca;
    ax.XAxis.Visible = 'on';
    yticks([5 10]);
    xticks(xintercepts);
    xticklabels(xlabels);
    xtickangle(ax, 20);
end

function A = apply_ylim(A)
    maxY = ones(size(A)) * 3;
    minY = ones(size(A)) * -1;
    A = bsxfun(@min, A, maxY);
    A = bsxfun(@max, A, minY);
end
