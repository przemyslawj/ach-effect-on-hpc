path = '/mnt/DATA/Clara/baseline/2018-09-06/signal';
path = '/mnt/DATA/Clara/diode_baseline/20190906/BD031ActiveStimON1_g0';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

animal_code = binName(1:2);
channelTable = findOmneticChannels(...
    '/mnt/DATA/Clara/diode_baseline/channels.csv',...
    animal_code);

secondOffset = 0;
meta = ReadMeta(binName, path);
keepChannels = intersect(meta.snsSaveChanSubset, channelTable.channel, 'sorted');
channelTable = channelTable(ismember(channelTable.channel, keepChannels), :);

lengthSeconds = min(120, str2double(meta.fileTimeSecs) - secondOffset);
nChans = meta.nChans;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = meta.nSamp / 2;
%dataArray = filter50Hz(dataArray, fs);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
dataArray = dataArray(ismember(meta.snsSaveChanSubset, keepChannels),:);
channel_std = std(dataArray, [], 2);

if isempty(channelTable)
  channelList = (1:nChans)';
end

showTraces(dataArray, fs, binName, channelTable);


