%%% Calculated std estimate for each mouse on the filtered signal from the whole day

datarootdir = '/mnt/DATA/chat_ripples/y-maze';

expstable.animal_date_group = strcat(...
    datestr(expstable.date, 'yyyy-mm-dd'), ...
    '_',...
    expstable.animal);
[g, gN] = grp2idx(expstable.animal_date_group); 

std_estimates = table();
for gi = 1:max(g)
    group_exps = find(g == gi);
    group_expstable = expstable(group_exps,:);
    
    signal_concat = [];
    trial_std = zeros(numel(group_exps), 2);
    for i = 1:numel(group_exps)

        animal_code = group_expstable.animal{i};
        state = group_expstable.state{i};
        dateddir = datestr(group_expstable.date(i), 'yyyy-mm-dd');
        signalpath = [ datarootdir filesep dateddir filesep strtrim(group_expstable.dirname{i})];
        binfile = dir([ signalpath filesep animal_code '*.bin']);
        if size(binfile, 1) == 0
            warning('No data files found at %s', signalpath);
            continue
        end
        fprintf('Processing file for date=%s file=%s\n', dateddir, binfile.name);
        meta = ReadMeta(binfile.name, binfile.folder);
        channels_file = reverse_channels_file;
        if ismember('reverse_channel_map', group_expstable.Properties.VariableNames) &&...
                (group_expstable.reverse_channel_map(i) == 0)
            channels_file = ord_channels_file;
        end
        channelTable = readChannelTable(...
            channels_file, animal_code, meta, selected_channels_only, use_diode);
        dataArray = ReadSGLXData(meta, secondOffset, str2double(meta.fileTimeSecs) - secondOffset);
        
        channelTable = channelTable(~strcmp(channelTable.location, 'Laser'), :);
        dataArray = dataArray(channelTable.rec_order,:);
        
        fs = 1250;
        dataArray = downsample(dataArray', round(meta.nSamp / fs))';
        emg_index = find(strcmp(channelTable.location, 'EMG'));
        keep_sample = excludeEMGNoisePeriods(dataArray(emg_index,:),...
                                             fs * 1);
        channelTable = channelTable(~strcmp(channelTable.location, 'EMG'), :);
        
        if use_diode
            [dataArray, channelTable] = subtractDiodeSignal(dataArray, channelTable);
        end
        
        trial_filtered = applyRippleFilter(dataArray, channelTable, fs);
        
        ripple_detection_signal = zeros(size(trial_filtered));
        for channel = 1:size(trial_filtered,1)
            ripple_detection_signal(channel,:) = GetRippleSignal(...
                (trial_filtered(channel, :))', fs);
        end

        ripple_detection_signal = ripple_detection_signal(:, keep_sample);
        std_estimate = std(ripple_detection_signal, [], 2)';
        %std_estimate = (median(ripple_detection_signal, 2) / 0.6745)';
        trial_std(i, :) = std_estimate;
        signal_concat = [signal_concat ripple_detection_signal];
    end
    
    nchannels = size(signal_concat, 1);
    std_estimate = median(signal_concat, 2) / 0.6745;
    %std_estimate = std(signal_concat, [], 2);
    channel_name = channelTable.channel_name;

    
    group_std_estimates = table(channel_name, std_estimate);
    group_std_estimates.date = repmat(group_expstable.date(1), nchannels, 1);
    group_std_estimates.animal = repmat(group_expstable.animal(1), nchannels, 1);
    group_std_estimates.animal_date_group = repmat(group_expstable.animal_date_group(1), nchannels, 1);
    
    std_estimates = [std_estimates; group_std_estimates];
end

%% Write to csv
filename_infix = '';
if use_diode
    filename_infix = '_diode';
end
if ~selected_channels_only
    filename_infix = [filename_infix '_all'];
end
outfile_suffix = [filename_infix '.csv'];
writetable(std_estimates, [datarootdir filesep 'std_estimates' outfile_suffix]);
