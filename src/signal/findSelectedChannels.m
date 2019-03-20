function [channelList] = findSelectedChannels(channelsFilepath, animal_code)

if isfile(channelsFilepath)
    electrodes = readtable(channelsFilepath);
    channelList = electrodes(strcmp(electrodes.animal, animal_code),:).channel;
else
    disp('No channels selected');
    channelList = [];
end

end
