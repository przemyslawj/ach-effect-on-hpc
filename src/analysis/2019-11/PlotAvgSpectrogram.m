datarootdir = '/mnt/DATA/chat_ripples/y-maze';
is_ymaze_trial = 1;
secondOffset = 2;
%datarootdir = '/mnt/DATA/chat_ripples/baseline';
%is_urethane_trial = 0;
%is_ymaze_trial = 0;
%secondOffset = 0;

animal_code = 'BS';
laserOn = 0;
date_str = '2019-11-13';
trials_fpath = [datarootdir filesep 'trials.csv'];
expstable = readtable(trials_fpath, 'ReadVariableNames', true);
expstable.dirname = strtrim(expstable.dirname);
reverse_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels_reversed.csv';
ord_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels.csv';
nexp = size(expstable, 1);
daytable = expstable(...
    strcmp(datestr(expstable.date, 'yyyy-mm-dd'), repmat({date_str}, nexp, 1)) &...
    strcmp(expstable.animal, repmat({animal_code}, nexp, 1)) &...
    expstable.laserOn == laserOn & expstable.correct == 1, :);
% daytable = expstable(...
%     strcmp(expstable.animal, repmat({animal_code}, nexp, 1)) &...
%     expstable.laserOn == laserOn & expstable.correct == 1, :);

result_table = table();

all_ripples = table();

all_bands = exp(0.7:0.05:5.55);

trial_period_lengths = [10 20 50 50];
%trial_period_lengths = [5 10 20 10];
% Matrix with spectrogram values indexed: channel location x 
% freq band index x common timestamp
band_zscore = zeros(2, numel(all_bands)-1, sum(trial_period_lengths));
% Count of trials: trial period x common timestamp
ntrials = zeros(1, sum(trial_period_lengths));
%daytable = daytable(1:2,:);

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
        else 
            continue;
        end
    end

    for chan_i = 1:size(channelTable, 1)
        loc = channelTable.location(chan_i);

        if strcmp(loc, 'EMG') || strcmp(loc, 'Laser')
            continue
        end
        chan_output_i = find(strcmp(loc, channelLocs));
        
        [wt, wfreqs] = cwt( dataArray(chan_i,:), 'morse', fs, 'ExtendSignal', true, ...
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
        for trial_period_i = 1:size(trialPeriods, 1)
            % Spectrogram signal
            start_index = max(1, trialPeriods.starts(trial_period_i));
            end_index = min(trialPeriods.ends(trial_period_i), size(dataArray,2));
            if end_index - start_index < 0
                continue
            end

            trial_period_offsets = cumsum([0 trial_period_lengths]);
            offset_i = trial_period_offsets(trial_period_i);
            period_len = trial_period_lengths(trial_period_i);
            ntrials_inc = zeros(1, size(ntrials, 2));
            period_pxx = band_power(:, start_index:end_index);
            % Compress the trial period length
            if size(period_pxx, 2) > period_len
                blockSize = [1, floor(size(period_pxx, 2) / period_len)];
                meanFilterFun = @(blockStruct) mean2(blockStruct.data(:));
                blockAvgSignal = blockproc(period_pxx, blockSize, meanFilterFun);
                blockAvgSignal = blockAvgSignal(:,1:period_len);
                ntrials_inc((offset_i + 1):(offset_i + period_len)) = 1;
            else % Keep the trial period length
                blockAvgSignal = zeros(size(period_pxx, 1), period_len);
                nsamples = size(period_pxx,2);
                blockAvgSignal(:,1:nsamples) = period_pxx;
                
                ntrials_inc((offset_i + 1):(offset_i + nsamples)) = 1;
            end
            ntrials = ntrials + ntrials_inc;
            periods_bandpower(:, (offset_i + 1):(offset_i + period_len)) = blockAvgSignal;
        end
        band_zscored_trial = z_score(periods_bandpower);
        %band_zscored_trial = movmean(band_zscored_trial, 3, 2);
        band_zscore(chan_output_i, :, :) = ...
            squeeze(band_zscore(chan_output_i, :, :)) + band_zscored_trial;
    end
end

%% Plot spectrograms
for chan_i = [1 2]
    f = figure('name', ['Channel=' channelLocs{chan_i} ],...
               'Position', [400 700 800 300]);
    zscore_mean = squeeze(band_zscore(chan_i, :, :)) ./ ...
        squeeze((ntrials / 2)); % Each trial counted twice: once per channel
    %zscore_mean = movmean(zscore_mean, 5, 2);
    zscore_mean = apply_ylim(zscore_mean);
    xintercepts = [1 cumsum(trial_period_lengths)];
    plotSpectrogram(zscore_mean, 1:sum(trial_period_lengths), ...
        all_bands(2:end), xintercepts);
    saveas(f, ['/home/prez/tmp/ymaze_chat_x_ai32/ymaze_spect' channelLocs{chan_i} '.svg'])
end

function [] = plotSpectrogram(wt_pow, time, wfreqs, xintercepts)
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
    xticklabels({'Start Zone','Maze Stem','Goal Zone', 'Last 10 sec'})
    xtickangle(ax, 20);
end

function A = apply_ylim(A)
    maxY = ones(size(A)) * 3;
    minY = ones(size(A)) * -1;
    A = bsxfun(@min, A, maxY);
    A = bsxfun(@max, A, minY);
end

function A = z_score(cfs)
    %cfs = log(cfs);
    A = cfs - mean(cfs, 2);
    A = bsxfun(@rdivide, A, std(A, [], 2));    
    %A = movmean(A, 10, 2);
end

function [] = draw_vlines(xintercepts, ylim)
    if isempty(xintercepts)
        return
    end
    for i = 1:numel(xintercepts)
        x = xintercepts(i);
        line([x, x], ylim, 'Color', 'white', 'LineWidth', 3);
    end
end
