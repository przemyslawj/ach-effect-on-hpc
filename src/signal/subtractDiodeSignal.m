function [data, diodeTable] = subtractDiodeSignal(...
        dataArray, channelTable)

data = [];
matchedChannels = [];
matchedChannelGroup = {};

matchedDataIndex = 1;
channelTable.channelGroup = strcat(channelTable.location, channelTable.side);
[g, gN] = grp2idx(channelTable.channelGroup); 
for gi = 1:max(g)
    loc_group_index = find(g == gi);
    if strcmp(channelTable.location{loc_group_index(1)}, 'EMG') || ...
            strcmp(channelTable.location{loc_group_index(1)}, 'Laser')
        data(matchedDataIndex,:) = dataArray(loc_group_index(1), :);
        channel = channelTable.channel(loc_group_index(1));
        channel_index = channelTable.channel_index(loc_group_index(1));
        matchedChannels(end + 1, :) = [channel, channel, channel_index, channel_index];
        matchedChannelGroup{end + 1} = channelTable.channelGroup{loc_group_index(1)};
        matchedDataIndex = matchedDataIndex + 1;
        continue
    end
    if numel(loc_group_index) < 2
        sprintf('Diode channel not found in %s', channelTable.location{loc_group_index})
        continue
    end
    
    locChannelTable = channelTable(loc_group_index,:);
    ind_comb = nchoosek(1:numel(loc_group_index), 2) ;
    for i = 1:size(ind_comb, 1)
        if locChannelTable.pair(ind_comb(i, 1)) == locChannelTable.pair(ind_comb(i, 2))
            continue
        end
        fst_channel = locChannelTable.channel(ind_comb(i, 1));
        snd_channel = locChannelTable.channel(ind_comb(i, 2));
        matchedChannelGroup{end + 1} = locChannelTable.channelGroup{1};
        fst_data_index = locChannelTable.channel_index(ind_comb(i, 1));
        snd_data_index = locChannelTable.channel_index(ind_comb(i, 2));
        matchedChannels(end + 1, :) = [fst_channel, snd_channel,...
            fst_data_index, snd_data_index];
        data(matchedDataIndex, :) = dataArray(fst_data_index, :) -...
                dataArray(snd_data_index, :);
        matchedDataIndex = matchedDataIndex + 1;
    end
    
end

diodeTable = array2table(matchedChannels, ...
    'VariableNames', {'channel1', 'channel2', 'channel1_index', 'channel2_index'});
diodeTable.channelGroup = matchedChannelGroup';
diodeTable.location = channelTable.location(matchedChannels(:,3));
ndiodes = size(diodeTable, 1);
diodeTable.channel_name = strcat(diodeTable.channelGroup, ...
    repmat(': ', ndiodes, 1),...
    num2str(diodeTable.channel1),...
    repmat('-', ndiodes, 1),...
    num2str(diodeTable.channel2));
diodeTable.channel = diodeTable.channel1;

end