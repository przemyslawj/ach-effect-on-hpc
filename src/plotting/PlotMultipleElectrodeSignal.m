%path = '/mnt/DATA/chat_ripples/baseline/';
path = '/mnt/DATA/chat_ripples/y-maze/';
%path = '/mnt/DATA/chat_ripples/y-maze/2019-11-13/signal';

[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

showFiltered = 1;
showDiode = 1;
reversed_channel_map = 1;

secondOffset = 0;
meta = ReadMeta(binName, path);

uppercase_idx = find(isstrprop(binName, 'upper'));
animal_code = binName(uppercase_idx(1:2));
    %'/mnt/DATA/chat_ripples/channel_desc/channels_reversed_baseline.csv',...
channelTable = readChannelTable(...
    '/mnt/DATA/chat_ripples/channel_desc/channels.csv',...
    animal_code, meta, reversed_channel_map);

lengthSeconds = min(120, str2double(meta.fileTimeSecs) - secondOffset);
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
nChans = meta.nChans;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = meta.nSamp / 10;
%dataArray = filter50Hz(dataArray, fs);
dataArray = dataArray(channelTable.rec_order,:);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
channel_std = std(dataArray, [], 2);

if isempty(channelTable)
  channelList = meta.snsSaveChanSubset;
end

if ~isempty(channelTable.location)
    for chan = 1:size(dataArray,1)    
        loc = channelTable.location{chan};
        if strcmp(loc, 'Laser')
            x = dataArray(chan,:);
            x = x - mean(x);
            dataArray(chan,:) = x / max(x) / 20;
        end
    end
end
shift = 5 * median(std(dataArray, [], 2));
showTraces(dataArray, fs, binName, channelTable, shift);

%% Filter LFP
filtered = applyRippleFilter(dataArray, channelTable, fs);

if showFiltered
    shift = 10 * median(std(filtered, [], 2));
    showTraces(filtered, fs, binName, channelTable, shift);
end

%% Diode signal
if showDiode
    [diodeData, diodeTable] = subtractDiodeSignal(dataArray, channelTable);
    showTraces(diodeData, fs, binName, diodeTable, shift);
end
if showDiode & showFiltered
    diodeFiltered = applyRippleFilter(diodeData, diodeTable, fs);
    showTraces(diodeFiltered, fs, binName, diodeTable, shift);
end
