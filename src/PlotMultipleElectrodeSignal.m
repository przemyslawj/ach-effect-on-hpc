path = 'E:\Ymaxe recording\2018-08-15';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

secondOffset = 0;
meta = ReadMeta(binName, path);

lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
fs = meta.nSamp;
nChans=21;
meta.nChans = nChans;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 1250;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

dataArray = filter50Hz(dataArray, fs);
channel_std = std(dataArray, [], 2);
%noisyChannels = find(channel_std > 0.1);
%dataArray(noisyChannels, :) = zeros(numel(noisyChannels), size(dataArray,2));

% shank 1L
channelList = 1:21;
showTraces(dataArray(channelList,:), fs, 'shank 1L', channelList);
