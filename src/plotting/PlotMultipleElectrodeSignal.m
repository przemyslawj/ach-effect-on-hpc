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
lengthSeconds = str2double(meta.fileTimeSecs) - secondOffset;
nChans = meta.nChans;

dataArray = ReadSGLXData(meta, secondOffset, lengthSeconds);

fs = meta.nSamp / 10;
%dataArray = filter50Hz(dataArray, fs);
dataArray = downsample(dataArray', round(meta.nSamp / fs))';
dataArray = dataArray(ismember(meta.snsSaveChanSubset, keepChannels),:);
channel_std = std(dataArray, [], 2);

if isempty(channelTable)
  channelList = (1:nChans)';
end

for chan = 1:size(dataArray,1)
    loc = channelTable.location{chan};
    if strcmp(loc, 'Laser')
        x = dataArray(chan,:);
        x = x - mean(x);
        dataArray(chan,:) = x / max(x) / 20;
    end
end
shift = 5 * median(std(dataArray, [], 2));
showTraces(dataArray, fs, binName, channelTable, shift);
%% Diode signal
[diodeData, diodeTable] = subtractDiodeSignal(dataArray, keepChannels, channelTable);
showTraces(diodeData, fs, binName, diodeTable, shift);

%% Filter LFP
filtered = zeros(size(dataArray));
passband = [80 250];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);

for chan = 1:size(dataArray,1)
    loc = channelTable.location{chan};
    x = dataArray(chan,:);
    x = x - mean(x);
    if strcmp(loc, 'EMG') 
        filtered(chan,:) = x / max(x) / 20;
    elseif strcmp(loc, 'Laser')
        filtered(chan,:) = x / 5;
    else
        filtered(chan,:) = filtfilt(b, a, x);
    end
end

shift = 10 * median(std(filtered, [], 2));
showTraces(filtered, fs, binName, channelTable, shift);

%% Diode signal
[filteredDiodeData, filteredDiodeTable] = subtractDiodeSignal(filtered, keepChannels, channelTable);
showTraces(filteredDiodeData, fs, binName, filteredDiodeTable, shift);
