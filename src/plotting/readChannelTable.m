function [channelTable] = readChannelTable(channels_csv, animal_code, ...
    meta, reversed_channel_map, selected_channels_only, use_diode)
%READCHANNELTABLE Returns table with channels and their description

if nargin < 5
    selected_channels_only = 0;
    use_diode = 0;
end
keep_emg = 1;

channelTable = findOmneticChannels(...
    channels_csv,...
    animal_code);
channelTable = channelTable(channelTable.reversed_map == reversed_channel_map, :);
keepChannels = intersect(meta.snsSaveChanSubset, channelTable.channel, 'sorted');
if isempty(keepChannels)
    keepChannels = (0:size(channelTable,1))';
end

channelTable = channelTable(ismember(channelTable.channel, keepChannels), :);
channelTable.rec_order = (find(ismember(meta.snsSaveChanSubset, keepChannels)))';

emg_index = find(strcmp(channelTable.location, 'EMG'), 1);
if selected_channels_only
    if use_diode
        selected_col = find(channelTable.is_pair_selected == 1);
    else
        selected_col = find(channelTable.is_channel_selected == 1);
    end
    if keep_emg && ~ismember(emg_index, selected_col)
        selected_col = [selected_col; emg_index];
    end
    channelTable = channelTable(selected_col,:);
end
channelTable.channel_index = (1:size(channelTable, 1))';

end

