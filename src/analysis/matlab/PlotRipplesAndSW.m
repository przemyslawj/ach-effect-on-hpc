max_recording_ripples = 100;

datarootdir = '/mnt/DATA/chat_ripples/baseline';
secondOffset = 0;
is_ymaze_trial = 0;
is_after_ymaze = 0;

diode_ripples = 1;
save_plots = 1;
selected_channels_only = 1;


channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels_reversed_ymaze.csv';
channels_file_gfp = '/mnt/DATA/chat_ripples/channel_desc/channels_gfp.csv';    
if ~is_ymaze_trial
    channels_file = '/mnt/DATA/chat_ripples/channel_desc/channels.csv';
    datarootdir = '/mnt/DATA/chat_ripples/baseline';
end

use_diode = 1;
ripples_filename = 'ripples_diode_th6.csv';

ripplestable = readtable([datarootdir filesep 'trial_results' filesep ripples_filename]);
if is_ymaze_trial
    secondOffset = 3;
    ripplestable = ripplestable(strcmp(ripplestable.stage_desc, 'DuringStim'),:);
end
ripplestable = ripplestable(strcmp(ripplestable.exp, 'main-effect' ) == 1, :);
%ripplestable = ripplestable(strcmp(ripplestable.animal, 'OS'),:);

trials_fpath = [datarootdir filesep 'trials.csv'];
if is_after_ymaze
    trials_fpath = [datarootdir filesep 'trials_after.csv'];
end
expstable = readtable(trials_fpath, 'ReadVariableNames', true);

ripplestable.dated_file_channel = strcat(...
    datestr(ripplestable.date, 'yyyy-mm-dd'), ...
    '_', ripplestable.file_name,...
    '_', ripplestable.channelLocation,...
    '_', 'laser', num2str(ripplestable.laserOn));
[g, gN] = grp2idx(ripplestable.dated_file_channel); 

for gi = 1:max(g)
    group_file_indecies = find(g == gi);
    group_ripplestable = ripplestable(group_file_indecies,:);    
    animal_code = group_ripplestable.animal{1};
    exp_name = group_ripplestable.exp{1};
    
    dateddir = datestr(group_ripplestable.date(1), 'yyyy-mm-dd');
    signalpath = [ datarootdir filesep dateddir filesep strtrim(group_ripplestable.file_name{1})];
    binfile = dir([ signalpath filesep '*.bin']);
    if size(binfile, 1) == 0
        warning('No data files found at %s', signalpath);
        continue
    end
    fprintf('Processing file for date=%s file=%s\n', dateddir, binfile.name);
    meta = ReadMeta(binfile.name, binfile.folder);
    trial_channels_file = channels_file;
    if strcmp(exp_name, 'main-effect') == 0
        trial_channels_file = channels_file_gfp;
    end
    expstable_row = find(expstable.date == group_ripplestable.date(1) & ...
      strcmp(expstable.animal, animal_code) & ...
      strcmp(expstable.dirname, group_ripplestable.file_name{1}) & ...
      strcmp(expstable.exp, exp_name));
    
    channelTable = readChannelTable(...
        trial_channels_file, animal_code, meta, ...
        expstable.reverse_channel_map(expstable_row),...
        selected_channels_only, use_diode);
    
    dataArray = ReadSGLXData(meta, secondOffset, str2double(meta.fileTimeSecs) - secondOffset);
    dataArray = dataArray(channelTable.rec_order,:);

    fs = 1250;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    time = (1:size(dataArray,2)) / fs;
    [dataDiode, channelTableDiode] = subtractDiodeSignal(dataArray, channelTable);
    filteredDiode = applyRippleFilter(dataDiode, channelTableDiode, fs);
    emgIdx = find(strcmp(channelTableDiode.location, 'EMG'));
    keep_sample_fewer = excludeEMGNoisePeriods(dataDiode(emgIdx,:), fs * 1);
    chan_name_parts = strrep(strsplit(group_ripplestable.channelName{1}, ' '), '-', '');
    channel_index_diode = find(...
        startsWith(strrep(channelTableDiode.channel_name, ' ', ''), ...
        strrep(group_ripplestable.channelName{1}, ' ', '')));
    if diode_ripples
        fst_chan = find(channelTable.channel == str2num(chan_name_parts{2}));
        snd_chan = find(channelTable.channel == str2num(chan_name_parts{3}));
    else
        fst_chan = find(strcmp(channelTable.channel_name,chan_name_parts{1}));
        snd_chan = fst_chan;
        channel_index_diode = fst_chan;
    end
    
    ripple_slots = 5;
    nslots_h = 1;
    nslots_w = ripple_slots / nslots_h;
    nripples = min(max_recording_ripples, size(group_ripplestable,1));
    if isempty(channel_index_diode)
        continue
    end
    
    batch_ripples = reshape(1:(ceil(nripples/ripple_slots) * ripple_slots), ...
                            [], ceil(nripples / ripple_slots));
    %% Plot single SWRs batch from recording
    for batch_i = 1:size(batch_ripples, 2)
        figName = sprintf('Channel %s, batch %d, laser %d', ...
            channelTable.channel_name{channel_index_diode}, ...
            batch_i,...
            group_ripplestable.laserOn(1));
        h = figure('name', figName,...
                   'Position', [100 1000 1200 500]);
        hAxFilteredDiode = [];
        hAxRawDiode = [];
        hAxRawChan1 = [];
        hAxRawChan2 = [];
        hAxEMG = [];
        for i = 1:ripple_slots
            ripple_i = batch_ripples(i, batch_i);
            if ripple_i > nripples
                break
            end
            channel_filtered = filteredDiode(channel_index_diode,:);
            channel_zscored_filtered = (channel_filtered - mean(channel_filtered)) / ...
                    std(channel_filtered(keep_sample_fewer));
            channel_signal_diode = dataDiode(channel_index_diode, :);

            peak_index = round(group_ripplestable.abs_peak_t(ripple_i) * fs) + 1;
            peak_offset_sec = group_ripplestable.abs_peak_t(ripple_i) - group_ripplestable.peak_t(ripple_i);
            halfwindow_len = int32(fs * 0.1);
            start_index = max(1, peak_index - halfwindow_len);
            end_index = min(peak_index + halfwindow_len, length(channel_filtered));
            indecies = start_index : end_index; 
            if isempty(indecies)
                continue
            end
            slot_h = floor((i - 1)/ nslots_w) + 1;
            slot_i = 2 * (slot_h - 1) * nslots_w + mod(i - 1, nslots_w) + 1;

            %subplot(2 * ripple_slots, slot_w, slot_h);
            hAxFilteredDiode(i) = subplot(nslots_h * 5, nslots_w, slot_i);
            %channel_zscored_filtered(indecies), fs, ...
            plotSWR(time(indecies), ...
                channel_filtered(indecies), fs, ...
                group_ripplestable.start_sec(i) + peak_offset_sec, ...
                group_ripplestable.end_time(i) + peak_offset_sec);

            hAxRawDiode(i) = subplot(nslots_h * 5, nslots_w, slot_i + nslots_w);
            plot(time(indecies), channel_signal_diode(indecies));
            %ylim([-0.3, 0.3]);

            hAxRawChan1(i) = subplot(nslots_h * 5, nslots_w, slot_i + 2 * nslots_w);
            plot(time(indecies), dataArray(fst_chan, indecies));
            hAxRawChan2(i) = subplot(nslots_h * 5, nslots_w, slot_i + 3 * nslots_w);
            plot(time(indecies), dataArray(snd_chan, indecies));
            hAxEMG(i) = subplot(nslots_h * 5, nslots_w, slot_i + 4 * nslots_w);
            plot(time(indecies), dataDiode(emgIdx, indecies));

        end
        if isempty(hAxRawDiode)
            continue
        end

        linkaxes(hAxRawDiode,'y');
        linkaxes(hAxFilteredDiode,'y');
        linkaxes(hAxRawChan1,'y');
        linkaxes(hAxRawChan2,'y');
        linkaxes(hAxEMG,'y');


        if save_plots
            fig_dir = fullfile(datarootdir, 'swrs', animal_code);
            if ~isfile(fig_dir)
                mkdir(fig_dir);
            end
            file_name = [ strrep(gN{gi}, '_signal/', '_') '-' num2str(batch_i) '.svg'];
            saveas(h, fullfile(fig_dir, file_name), 'svg');
        else
            while true
            w = waitforbuttonpress; 
            switch w 
                case 1 % (keyboard press) 
                  key = get(gcf,'currentcharacter'); 
                      switch key
                          case 110 % n key
                              break
                          case 27 % escape key
                              return
                          otherwise 
                              % Wait for a different command. 
                      end
                end
            end
        end
        close(h)
    end
end