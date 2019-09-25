path = '/mnt/DATA/Clara/baseline/2018-09-06/signal';
path = '/mnt/DATA/Clara/diode_baseline/20190906/BD031ActiveStimON1_g0';
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
% 
% % shank 1L
% channelList = [6 5 9 12 7 2 0 1 15] + 1;
% showTraces(dataArray, fs, binName, channelList');
% 
% % shank 1R
% channelList = [10 11 4 3 13 14 8 15] + 1;
% showTraces(dataArray, fs, binName, channelList');
% 
% % shank 2L
% channelList = [25 26 22 19 24 29 31 30 16] + 1;
% showTraces(dataArray, fs, binName, channelList');
% 
% % shank 2R
% channelList = [21 20 27 28 18 17 23 16] + 1;
% showTraces(dataArray, fs, binName, channelList');

