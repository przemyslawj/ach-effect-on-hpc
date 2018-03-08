%% setup results
datarootdir = '/media/prez/DATA/Prez/N&A_rest';
statefilepath = [datarootdir filesep 'state.csv'];
statetable = readtable(statefilepath, 'ReadVariableNames', true);

nexp = size(statetable, 1);
psd_table = statetable;
nChans = 32;
psd_table.pow_slow = zeros(nexp, nChans);
psd_table.pow_theta = zeros(nexp, nChans);
psd_table.pow_slow_gamma = zeros(nexp, nChans);
psd_table.pow_above_theta = zeros(nexp, nChans);
psd_table.has_ripples = zeros(nexp, nChans);
psd_table.swr_freq = zeros(nexp, nChans);
psd_table.swr_powerhz = zeros(nexp, nChans);
psd_table.ripple_peakpow = zeros(nexp, nChans);
psd_table.ripple_dur = zeros(nexp, nChans);
psd_table.dom_freq = zeros(nexp, nChans);

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
    sec_length = sec_end - sec_start;
    
    % calculate PSD only for trial_duration
    psd_sec_start = 0;
    psd_sec_end = sec_end - sec_start;
    trial_duration = 10;
    if strcmp(statetable.laser{i}, 'OFF')
        if sec_start == 0
            psd_sec_start = sec_end - sec_start - trial_duration;
        else
            psd_sec_end = min(sec_length, psd_sec_start + trial_duration);
        end
    else
        psd_sec_end = min(sec_length, psd_sec_start + trial_duration);
    end
    
    dataArray = ReadSGLXData(meta, sec_start, sec_length);

    fs = 1250;
    time = (1:size(dataArray,2)) / fs;
    dataArrayOrg = dataArray;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    nChans = size(dataArray, 1);

    dataArray = filter50Hz(dataArray, fs);
    channel_std = std(dataArray, [], 2);
    downfs = 125;
    downsampledDataArray = downsample(dataArray', round(fs / downfs))';
    
    % Filter LFP for SWR detection
    filtered = zeros(size(dataArray));
    passband = [100 250];
    nyquist = fs / 2;
    filterOrder = 4;
    filterRipple = 20;
    [b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);
    for chan = 1:nChans
        filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
    end

    %% calculate PSD

    slow = [0 4];
    theta = [5 10];
    slow_gamma = [35 45];
    above_theta = [11, (fs / 2 - 1)];

    for channel = 1:nChans
        if channel_std(channel) > 0.1
            continue
        end
        
        %[pxx, freqs] = pwelch(dataArray(channel,:), ...
        %    floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
        psd_x = downsampledDataArray(channel,:);
        psd_xx = psd_x((floor(downfs * psd_sec_start) + 1) : min(floor(downfs * psd_sec_end), size(psd_x,2)));
        if isempty(psd_xx)
            continue
        end
        
        [cws, freqs] = cwt(psd_xx, 'morse', downfs);
        
        pxx = mean(abs(cws .^ 2), 2);
        [maxVal, maxValIndex] = max(pxx);
        psd_table.dom_freq(i, channel) = freqs(maxValIndex);
        
        psd_table.pow_slow(i, channel) = TotalBandPower2(freqs, pxx, slow);
        psd_table.pow_theta(i, channel) = TotalBandPower2(freqs, pxx, theta);
        psd_table.pow_slow_gamma(i, channel) = TotalBandPower2(freqs, pxx, slow_gamma);
        psd_table.pow_above_theta(i, channel) = TotalBandPower2(freqs, pxx, above_theta);
        

        [ripples, sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(channel,:)', ...
                             'frequency', fs, ...
                             'thresholds', [2 4 0.02],...
                             'durations', [30 20 300]);
         
         if ~isempty(ripples)
            psd_table.has_ripples(i, channel) = any(ripples(:,2) > psd_sec_start);
            psd_table.swr_freq(i, channel) = size(ripples,1) / sec_length;
            psd_table.ripple_peakpow(i, channel) = mean(ripples(:,4));
            psd_table.ripple_dur(i, channel) = mean(ripples(:,3) - ripples(:,1));
         end
         psd_table.swr_powerhz(i, channel) = ...
             mean(abs(normalizedSquaredSignal));
    end
end

writetable(psd_table, 'psd_table2.csv');
            
function [ pow ] = TotalBandPower(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,band,'psd');
end

function [ pow ] = TotalBandPower2(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    band_freqs = find(f1 >= band(1) & f1 <= band(2));
    pow = -trapz(f1(band_freqs), log10(pxx(band_freqs))) / (band(2) - band(1));
end