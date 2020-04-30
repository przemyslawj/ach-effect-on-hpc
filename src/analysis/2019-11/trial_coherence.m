%% setup results
use_diode = 1;
selected_channels_only = 1;
is_urethane = 0;
is_after_ymaze = 0;
%datarootdir = '/mnt/DATA/chat_ripples/y-maze';
%is_ymaze_trial = 1;
%secondOffset = 2.5;
datarootdir = '/mnt/DATA/chat_ripples/baseline';
is_ymaze_trial = 0;
secondOffset = 0;

trials_fpath = [datarootdir filesep 'trials.csv'];
expstable = readtable(trials_fpath, 'ReadVariableNames', true);
expstable.dirname = strtrim(expstable.dirname);
reverse_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels_reversed_same_side.csv';
ord_channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels.csv';

nexp = size(expstable, 1);
result_table = table();

% Predefined frequencies
all_bands = exp(0.1:0.035:5.6);
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

    fs = 5000;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    time = (1:size(dataArray,2)) / fs;
    if use_diode
        [dataArray, channelTable] = subtractDiodeSignal(dataArray, channelTable);
    end
    
    nchans = size(dataArray, 1);

    laserChannelIdx = find(strcmp(channelTable.location, 'Laser'));
    emgIdx = find(strcmp(channelTable.location, 'EMG'));
    if is_ymaze_trial == 0
        trialPeriods = extractTrialPeriodsFromLaser(dataArray,...
            laserChannelIdx, fs, 100);
        
        trialPeriods = trialPeriods(~strcmp(trialPeriods.stage_desc, 'after_stim'), :);
        % assumes trials ordered
        trialPeriods.trial_ordinal = (repelem(1:(size(trialPeriods,1)/2), 1, 2))';
    else
        tracking_filepath = [ datarootdir filesep get_trackingfilepath(dateddir, binfile.name)];
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
            trialPeriods = [];
        end
    end
    
    if is_after_ymaze
        trialPeriods = table();
        trialPeriods.starts = 0;
        trialPeriods.ends = int32(str2double(meta.fileTimeSecs) * fs);
        trialPeriods.laserOn = expstable.laserOn(i);
        trialPeriods.stage_desc = {'after'};
    end
    
    if isempty(trialPeriods)
        warning('No trial periods found in the recording')
        continue
    end
    

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
            theta = [5 10];
            above_theta = [11 20];
            slow_gamma = [25 45];
        end
    end
    
    ca_channels = find(startsWith(channelTable.location, 'CA'));
    
    signal_pair = dataArray(ca_channels, :);
    
    for trial_ordinal = (unique(trialPeriods.trial_ordinal))'
        trialSections = trialPeriods(trialPeriods.trial_ordinal == trial_ordinal, :);
        
        trial_start = max(1, min(trialSections.starts));
        trial_end = min(max(trialSections.ends), size(dataArray,2));
        trial_signal = signal_pair(:, trial_start : trial_end);
        keep_trial_sample = intersect(...
            excludeEMGNoisePeriods(trial_signal(ca_channels(1),:), fs * 0.2, 3),...
            excludeEMGNoisePeriods(trial_signal(ca_channels(2),:), fs * 0.2, 3));

        sec_start = max(0, double(trial_start) / fs);
        sec_end = double(trial_end) / fs;
        sec_length = sec_end - sec_start;

        if sec_length <= min_section_dur_sec
            continue
        end

        %[wcoh, wcs, freqs] = wcoherence(trial_signal(1,:), trial_signal(2,:), fs);
        %freqs = fliplr(freqs')';
        %wcoh = fliplr(wcoh')';
        
        for section_i = 1:size(trialSections,1)
            section_start_i = max(1, trialSections.starts(section_i) - trial_start + 1);
            section_end_i = min(trialSections.ends(section_i) - trial_start + 1,...
                               size(trial_signal,2));
            section_indecies = section_start_i:section_end_i;
            keep_section_samples = intersect(section_indecies, keep_trial_sample);

            if numel(keep_section_samples) < numel(section_indecies)
                nrejected_indecies = numel(section_indecies) - numel(keep_section_samples);
                warning(['Rejected ' num2str(nrejected_indecies / fs) ' s out of '...
                    num2str(numel(section_indecies) / fs) ' s due to noise'])
            end

            if numel(keep_section_samples) < 0.2 * fs
                warning(['Section ' num2str(section_i) ' too short, skippnig'])
                continue
            end                           
                          
            [section_coh, freqs] = mscohere(...
                trial_signal(1,section_indecies),...
                trial_signal(2,section_indecies), ...
                hanning(fs/5), ceil(fs/10), all_bands, fs);
            
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            result_table.date(entry_i) = {dateddir};
            result_table.animal(entry_i) = {animal_code};
            result_table.file_name(entry_i) = expstable.dirname(i);
            trial = [expstable.dirname{i} '_' num2str(trial_ordinal)];
            result_table.trial(entry_i) = {trial};
            result_table.stage_desc(entry_i) = trialSections.stage_desc(section_i);
            result_table.laserOn(entry_i) = trialSections.laserOn(section_i);
            
            result_table.coh_slow(entry_i) = CoherenceMean2(freqs, section_coh, slow);
            result_table.coh_theta(entry_i) = CoherenceMean2(freqs, section_coh, theta);
            result_table.coh_slow_gamma(entry_i) = CoherenceMean2(freqs, section_coh, slow_gamma);
            result_table.coh_med_gamma(entry_i) = CoherenceMean2(freqs, section_coh, med_gamma);
            result_table.coh_fast_gamma(entry_i) = CoherenceMean2(freqs, section_coh, fast_gamma);
            result_table.coh_above_theta(entry_i) = CoherenceMean2(freqs, section_coh, above_theta);
            result_table.coh_ripple_band(entry_i) = CoherenceMean2(freqs, section_coh, ripple_band);
            
            for j=1:(numel(all_bands)-1)
                result_table.all_coh(entry_i, j) = section_coh(j);
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

outfile_suffix = [filename_infix '.csv'];
writetable(result_table, [datarootdir filesep 'coherence_table' outfile_suffix]);

function [ mcoh ] = CoherenceMean(freqs, wcoh, band)
    band_freqs = freqs >= band(1) & freqs <= band(2);
    mcoh = mean(mean(wcoh(band_freqs,:),2));
end

function [ pow ] = CoherenceMean2(freqs, wcoh, band)
    band_freqs_i = find(freqs >= band(1) & freqs <= band(2));
    band_freqs = freqs(band_freqs_i);
    pow = trapz(band_freqs, wcoh(band_freqs_i)) / (max(band_freqs) - min(band_freqs));
end

