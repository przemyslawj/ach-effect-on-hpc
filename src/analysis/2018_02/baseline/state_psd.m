%% setup results
datarootdir = '/media/prez/DATA/Prez/N&A_rest';
%datarootdir = '/media/prez/DATA/Prez/N&A_rest/2018-03-13/signal';
statefilepath = [datarootdir filesep 'state.csv'];
statetable = readtable(statefilepath, 'ReadVariableNames', true);

nexp = size(statetable, 1);
psd_table = statetable;
nChans = 32;
psd_table.pow_slow = zeros(nexp, nChans);
psd_table.pow_theta = zeros(nexp, nChans);
psd_table.pow_slow_gamma = zeros(nexp, nChans);
psd_table.pow_med_gamma = zeros(nexp, nChans);
psd_table.pow_fast_gamma = zeros(nexp, nChans);
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
    
    if sec_length == 0
        continue
    end
    
    trial_duration = 10;
    %trial_duration = sec_length;
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
    downfs = 625;
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

    slow = [0.8 5];
    theta = [6 10];
    slow_gamma = [30 80];
    med_gamma = [60 120];
    fast_gamma = [100 250];
    above_theta = [10.5, 250];

    for channel = 1:nChans
        if channel_std(channel) > 0.1
            continue
        end
        
        psd_x = downsampledDataArray(channel,:);
        psd_xx = psd_x((floor(downfs * psd_sec_start) + 1) : min(floor(downfs * psd_sec_end), size(psd_x,2)));
        if isempty(psd_xx)
            continue
        end
        
        %[pxx, freqs] = pwelch(psd_xx, ...
        %    floor(downfs / 4), floor(downfs / 8), floor(downfs / 2), downfs);
        [cws, freqs] = cwt(psd_xx, 'morse', downfs);        
        pxx = mean(abs(cws .^ 2), 2);
        freqs = fliplr(freqs')';
        pxx = fliplr(pxx')';
        
        [maxVal, maxValIndex] = max(pxx);
        psd_table.dom_freq(i, channel) = freqs(maxValIndex);
        
        psd_table.pow_slow(i, channel) = TotalBandPower(freqs, pxx, slow);
        psd_table.pow_theta(i, channel) = TotalBandPower(freqs, pxx, theta);
        psd_table.pow_slow_gamma(i, channel) = TotalBandPower(freqs, pxx, slow_gamma);
        psd_table.pow_med_gamma(i, channel) = TotalBandPower(freqs, pxx, med_gamma);
        psd_table.pow_fast_gamma(i, channel) = TotalBandPower(freqs, pxx, fast_gamma);
        psd_table.pow_above_theta(i, channel) = TotalBandPower(freqs, pxx, above_theta);
        

        [ripples, sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(channel,:)', ...
                             'frequency', fs, ...
                             'thresholds', [2 4.5 0.01],...
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

writetable(psd_table, 'psd_table_rest_cwt2.csv');
            
function [ pow ] = TotalBandPower(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,band,'psd');
end

function [ pow ] = TotalBandPower2(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    band_freqs = find(f1 >= band(1) & f1 <= band(2));
    pow = -trapz(f1(band_freqs), log10(pxx(band_freqs))) / (band(2) - band(1));
end