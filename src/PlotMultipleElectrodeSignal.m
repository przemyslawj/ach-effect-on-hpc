path = '/media/prez/DATA/Prez/N&A_rest/';
[binName, path] = uigetfile('*.bin', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

secondOffset = 0;
meta = ReadMeta(binName, path);

lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
fs = meta.nSamp;
nChans=32;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = 1250;
dataArray = downsample(dataArray', round(meta.nSamp / fs))';

dataArray = filter50Hz(dataArray, fs);
channel_std = std(dataArray, [], 2);
noisyChannels = find(channel_std > 0.1);
dataArray(noisyChannels, :) = zeros(numel(noisyChannels), size(dataArray,2));

% shank 1L
channelList = [6 5 9 12 7 2 0 1 15] + 1;
showTraces(dataArray(channelList,:), fs, 'shank 1L', channelList);

% shank 1R
channelList = [10 11 4 3 13 14 8 15] + 1;
showTraces(dataArray(channelList,:), fs, 'shank 1R', channelList);

% shank 2L
channelList = [25 26 22 19 24 29 31 30 16] + 1;
showTraces(dataArray(channelList,:), fs, 'shank 2L', channelList);

% shank 2R
channelList = [21 20 27 28 18 17 23 16] + 1;
showTraces(dataArray(channelList,:), fs, 'shank 2R', channelList);
