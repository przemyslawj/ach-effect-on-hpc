function [channelTable] = findOmneticChannels(channelsFilepath, animal_code)

if isfile(channelsFilepath)
    electrodes = readtable(channelsFilepath);
    channelTable = electrodes(strcmp(electrodes.animal_code, animal_code),:);
    channelTable = sortrows(channelTable, {'channel'});
    channelTable.channel_index = (1:size(channelTable, 1))';
    channelTable.channel_name = strcat(channelTable.location, ...
                                       channelTable.side, ...
                                       num2str(channelTable.pair));
else
    warning(['File not found for channels at ' channelsFilepath]);
    channelTable = table();
end

end
