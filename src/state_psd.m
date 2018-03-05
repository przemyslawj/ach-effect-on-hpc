%% setup results
datarootdir = '/media/prez/DATA/Prez/N&A_rest/';
statefilepath = [datarootdir filesep 'state.csv'];
statetable = readtable(statefilepath, 'ReadVariableNames', true);

nexp = size(statetable, 1);
psd_table = statetable;
nChans = 32;
psd_table.pow_slow = zeros(nexp, nChans);
psd_table.pow_theta = zeros(nexp, nChans);
psd_table.pow_slow_gamma = zeros(nexp, nChans);
psd_table.pow_above_theta = zeros(nexp, nChans);

%% process single experiments
for i = 1:nexp
    binfile = dir([datarootdir filesep statetable.experiment{i} '*.bin']);
    fprintf('Processing file: %s, experiment: %s \n', ...
        binfile.name, statetable.experiment{i});
    meta = ReadMeta(binfile.name, binfile.folder);

    sec_start = statetable.time_sec(i);
    sec_end = str2double(meta.fileTimeSecs);
    if i + 1 < size(statetable, 1) ...
            && strcmp(statetable.experiment{i}, statetable.experiment{i+1})
        sec_end = statetable.time_sec(i+1);
    end
    % skip first 3sec after change of the laser state
    if sec_start > 0
        sec_start = sec_start + 3;
    end

    dataArray = ReadSGLXData(meta, sec_start, sec_end - sec_start);

    fs = 1250;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    nChans = size(dataArray, 1);

    dataArray = filter50Hz(dataArray, fs);
    channel_std = std(dataArray, [], 2);

    %% calculate PSD

    slow = [0 4];
    theta = [5 10];
    slow_gamma = [35 45];
    above_theta = [11, (fs / 2 - 1)];


    for channel = 1:nChans
        if channel_std(channel) > 0.1
            continue
        end
        
        [pxx, freqs] = pwelch(dataArray(channel,:), ...
            floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
        psd_table.pow_slow(i, channel) = TotalBandPower(freqs, pxx, slow);
        psd_table.pow_theta(i, channel) = TotalBandPower(freqs, pxx, theta);
        psd_table.pow_slow_gamma(i, channel) = TotalBandPower(freqs, pxx, slow_gamma);
        psd_table.pow_above_theta(i, channel) = TotalBandPower(freqs, pxx, above_theta);
    end
end

writetable(psd_table, 'psd_table.csv');
save('psd_table.mat', 'psd_table');
            
function [ pow ] = TotalBandPower(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,band,'psd');
end