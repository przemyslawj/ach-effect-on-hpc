path = '/mnt/DATA/Clara/baseline/2018-09-06/signal';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

animal_code = binName(1:2);
channelList = findSelectedChannels(...
    '/mnt/DATA/Clara/ymaze/selected_electrodes.csv',...
    animal_code);

secondOffset = 0;
meta = ReadMeta(binName, path);

lengthSeconds = min(120, str2double(meta.fileTimeSecs) - secondOffset);
fs = meta.nSamp;
nChans = meta.nChans;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 2 * 1250;
fs = meta.nSamp / 2;
%dataArray = filter50Hz(dataArray, fs);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
channel_std = std(dataArray, [], 2);

if isempty(channelList)
  channelList = (1:nChans)';
end
showTraces(dataArray, fs, binName, channelList);

