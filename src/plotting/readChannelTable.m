function [channelTable] = readChannelTable(channels_csv, animal_code, meta)
%READCHANNELTABLE Returns table with channels and their description

channelTable = findOmneticChannels(...
    channels_csv,...
    animal_code);
keepChannels = intersect(meta.snsSaveChanSubset, channelTable.channel, 'sorted');
if isempty(keepChannels)
    keepChannels = (0:nChans)';
end
channelTable = channelTable(ismember(channelTable.channel, keepChannels), :);
channelTable.rec_order = (find(ismember(meta.snsSaveChanSubset, keepChannels)))';
end

