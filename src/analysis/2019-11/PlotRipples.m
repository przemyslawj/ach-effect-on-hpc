datarootdir = '/mnt/DATA/chat_ripples/y-maze';
secondOffset = 3;

save_plots = 1;

ripples_filename = 'ripples_diode_th6.csv';
ripplestable = readtable([datarootdir filesep ripples_filename]);
ripplestable = ripplestable(strcmp(ripplestable.stage_desc, 'DuringStim'),:);

ripplestable.dated_file_channel = strcat(...
    datestr(ripplestable.date, 'yyyy-mm-dd'), ...
    '_', ripplestable.file_name,...
    '_', ripplestable.channelLocation);
[g, gN] = grp2idx(ripplestable.dated_file_channel); 

for gi = 1:max(g)
    group_file_indecies = find(g == gi);
    group_ripplestable = ripplestable(group_file_indecies,:);    
    animal_code = group_ripplestable.animal{1};
    
    dateddir = datestr(group_ripplestable.date(1), 'yyyy-mm-dd');
    signalpath = [ datarootdir filesep dateddir filesep strtrim(group_ripplestable.file_name{1})];
    binfile = dir([ signalpath filesep animal_code '*.bin']);
    if size(binfile, 1) == 0
        warning('No data files found at %s', signalpath);
        continue
    end
    fprintf('Processing file for date=%s file=%s\n', dateddir, binfile.name);
    meta = ReadMeta(binfile.name, binfile.folder);
    channels_file = reverse_channels_file;
    if ismember('reverse_channel_map', group_ripplestable.Properties.VariableNames) &&...
            (expstable.reverse_channel_map(1) == 0)
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
    filtered = applyRippleFilter(dataArray, channelTable, fs);
    emgIdx = find(strcmp(channelTable.location, 'EMG'));
    keep_sample_fewer = excludeEMGNoisePeriods(dataArray(emgIdx,:), fs * 1);
    
    channel_index = find(...
        strcmp(strrep(channelTable.channel_name, ' ', ''), ...
        strrep(group_ripplestable.channelName{1}, ' ', '')));
    figName = sprintf('Channel %s', channelTable.channel_name{channel_index});
    h = figure('name', figName,...
               'Position', [100 1000 1200 500]);
    ripple_slots = 8;
    nslots_h = 2;
    nslots_w = ripple_slots / nslots_h;
    nripples = min(ripple_slots, size(group_ripplestable,1));
    hAxFiltered = [];
    hAxRaw = [];
    for i=1:nripples
        channel_filtered = filtered(channel_index,:);
        channel_zscored_filtered = (channel_filtered - mean(channel_filtered)) / ...
                std(channel_filtered(keep_sample_fewer));
        channel_signal = dataArray(channel_index, :);
        
        peak_index = round(group_ripplestable.abs_peak_t(i) * fs) + 1;
        peak_offset_sec = group_ripplestable.abs_peak_t(i) - group_ripplestable.peak_t(i);
        halfwindow_len = int32(fs * 0.1);
        start_index = max(1, peak_index - halfwindow_len);
        end_index = min(peak_index + halfwindow_len, length(channel_filtered));
        indecies = start_index : end_index; 
        slot_h = floor((i - 1)/ nslots_w) + 1;
        slot_i = 2 * (slot_h - 1) * nslots_w + mod(i - 1, nslots_w) + 1;
        
        %subplot(2 * ripple_slots, slot_w, slot_h);
        hAxFiltered(i) = subplot(nslots_h * 2, nslots_w, slot_i);
        plotSWR(time(indecies), channel_zscored_filtered(indecies), fs, ...
            group_ripplestable.start_sec(i) + peak_offset_sec, ...
            group_ripplestable.end_time(i) + peak_offset_sec);

        ylim_val = 8;
        %ylim([-ylim_val ylim_val]);

        hAxRaw(i) = subplot(nslots_h * 2, nslots_w, slot_i + nslots_w);
        %subplot(2 * ripple_slots, slot_w, slot_h + 1);
        plot(time(indecies), channel_signal(indecies));
        %ylim([-0.3, 0.3]);
        
    end
    linkaxes(hAxRaw,'y');
    linkaxes(hAxFiltered,'y');
    
    if save_plots
        fig_dir = fullfile(datarootdir, 'swrs100_250', animal_code);
        if ~isfile(fig_dir)
            mkdir(fig_dir);
        end
        file_name = [ strrep(gN{gi}, '_signal/', '_') '.svg'];
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