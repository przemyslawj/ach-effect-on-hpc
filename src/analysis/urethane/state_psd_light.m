%% setup results
rootdir = '/mnt/DATA/Clara/urethane/';
experiment = 'inhibition';
datarootdir = [ rootdir experiment filesep ];

txtfiles = dir([datarootdir filesep '*.txt']);

%expstable = expstable(strcmp(expstable.animal, animal_code),:);

nchans = 1;

nexp = size(txtfiles, 1);
result_table = table();

all_ripples = table();

result_table.animal = repmat('XXXXXX', nexp * 3 * nchans, 1);

result_table.row = zeros(nexp * 3 * nchans, 1);
result_table.trial = zeros(nexp * 3 * nchans, 1);
result_table.laserOn = zeros(nexp * 3 * nchans, 1);
result_table.channel = zeros(nexp * 3 * nchans, 1);
all_bands = [0.1 0.5 1 1.5, exp(0.7:0.05:5.3)];
result_table.all_psd_xx = zeros(nexp * 3 * nchans, numel(all_bands));

result_table.pow_slow = zeros(nexp * 3 * nchans, 1);
result_table.pow_theta = zeros(nexp * 3 * nchans, 1);
result_table.peak_theta = zeros(nexp * 3 * nchans, 1);

result_table.pow_slow_gamma = zeros(nexp * 3  * nchans, 1);
result_table.pow_med_gamma = zeros(nexp * 3 * nchans, 1);
result_table.pow_fast_gamma = zeros(nexp * 3 * nchans, 1);
result_table.pow_above_theta = zeros(nexp * 3 * nchans, 1);
result_table.has_ripples = zeros(nexp * 3 * nchans, 1);
result_table.swr_incidence = zeros(nexp * 3 * nchans, 1);
result_table.swr_powerhz = zeros(nexp * 3 * nchans, 1);
result_table.ripple_peakpow = zeros(nexp * 3 * nchans, 1);
result_table.ripple_dur = zeros(nexp * 3 * nchans, 1);
result_table.dom_freq = zeros(nexp * 3 * nchans, 1);

%% process single experiments
entry_i = 1;
for i = 1:nexp

    fs = 1000;
    nChans = 1;
    channelList = 1;

    fname = [datarootdir filesep txtfiles(i).name]
    animal = txtfiles(i).name(1:6);
    datatable = readtable(fname);
    dataArray = datatable{:,1};
    dataArray = dataArray';
    lengthSeconds = size(dataArray, 2) / fs;

    rec_times = [0 15;
                 15 30;
                 45 60; ];    
%     rec_times = [0 60;
%                  60 120;
%                  120 180; ];

    time = (1:size(dataArray,2)) / fs;

    %dataArray = filter50Hz(dataArray, fs);
    channel_std = std(dataArray, [], 2);
    downfs = 500;
    downsampledDataArray = downsample(dataArray', round(fs / downfs))';

    % Filter LFP for SWR detection
    filtered = zeros(size(dataArray));
    passband = [80 250];
    nyquist = fs / 2;
    filterOrder = 4;
    filterRipple = 20;
    [b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);
    for chan = 1:nchans
        filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
    end

    %% calculate PSD

    slow = [0.1 2];
    theta = [2 5];
    slow_gamma = [20 40];
    med_gamma = [60 80];
    fast_gamma = [80 150];
    above_theta = [7, 20];

    for channel = 1:nchans
        sd = 0;      
        for rec_times_i = 1:size(rec_times, 1)
            sec_start = rec_times(rec_times_i, 1);
            sec_end = rec_times(rec_times_i, 2);
            sec_length = sec_end - sec_start;

            if sec_length <= 0
                entry_i = entry_i + 1;
                continue
            end
            
            result_table.row(entry_i) = entry_i;
            result_table.laserOn(entry_i) = mod(rec_times_i - 1, 2);
            result_table.animal(entry_i,:) = animal;
            result_table.trial(entry_i) = i;
            result_table.laserOn(entry_i) = mod(rec_times_i - 1, 2);

            psd_xx = downsampledDataArray(channel,sec_start * downfs + 1:sec_end * downfs);

            %[pxx, freqs] = pwelch(psd_xx, ...
            %    floor(downfs / 4), floor(downfs / 8), floor(downfs / 2), downfs);
            [cws, freqs] = cwt(psd_xx, 'amor', downfs);        
            pxx = median(abs(cws .^ 2), 2);
            freqs = fliplr(freqs')';
            pxx = fliplr(pxx')';

            [maxVal, maxValIndex] = max(pxx);
            result_table.dom_freq(entry_i) = freqs(maxValIndex);

            result_table.pow_slow(entry_i) = TotalBandPower(freqs, pxx, slow);
            result_table.pow_theta(entry_i) = TotalBandPower(freqs, pxx, theta);
            result_table.peak_theta(entry_i) = PeakFreq(freqs, pxx, theta);
            result_table.pow_slow_gamma(entry_i) = TotalBandPower(freqs, pxx, slow_gamma);
            result_table.pow_med_gamma(entry_i) = TotalBandPower(freqs, pxx, med_gamma);
            result_table.pow_fast_gamma(entry_i) = TotalBandPower(freqs, pxx, fast_gamma);
            result_table.pow_above_theta(entry_i) = TotalBandPower(freqs, pxx, above_theta);
            for j=1:(numel(all_bands)-1)
                result_table.all_psd_xx(entry_i, j) = TotalBandPower(freqs, pxx, [all_bands(j) all_bands(j+1)]);
            end

             filtered_start_i = sec_start * fs + 1;
             filtered_end_i = min(sec_end * fs + 1, size(filtered,2));
             [section_ripples, sd, normalizedSquaredSignal] = MyFindRipples(...
                     (time(filtered_start_i : filtered_end_i))', ...
                     (filtered(channel, filtered_start_i : filtered_end_i))', ...
                     'frequency', fs, ...
                     'thresholds', [2 4.0 0.01], ...
                     'durations', [10 30 400], ....
                     'std', sd);
             nripples = size(section_ripples, 1);
             if ~isempty(section_ripples)
                result_table.has_ripples(entry_i) = 1;
                result_table.swr_incidence(entry_i) = size(section_ripples,1) / sec_length;
                result_table.ripple_peakpow(entry_i) = mean(section_ripples(:,4));
                result_table.ripple_dur(entry_i) = mean(section_ripples(:,3) - section_ripples(:,1));
                result_table.ripple_freq(entry_i) = mean(section_ripples(:,5));
               
                ripples_table = array2table(section_ripples);
                ripples_table.Properties.VariableNames = {'start_sec', 'peak_t', 'end_time', 'peakpow', 'peak_freq'};
                ripples_table.trial = repmat(i, nripples, 1);
                ripples_table.laserOn = repmat(result_table.laserOn(entry_i), nripples, 1);
 
                all_ripples = [all_ripples; ripples_table];
             end
             result_table.swr_powerhz(entry_i) = mean(abs(normalizedSquaredSignal));
            
             entry_i = entry_i + 1;
        end
       
    end
end

writetable(result_table, [rootdir filesep 'psd_table_urethane_' experiment '.csv']);
writetable(all_ripples, [rootdir filesep 'ripples_urethane_' experiment '.csv']);
            
function [ pow ] = TotalBandPower(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    pow = bandpower(pxx,f1,band,'psd');
end

function [ pow ] = TotalBandPower2(f1, pxx, band)
    %pow = sum(10 * log10(pxx(f1 >= band(1) & f1 <= band(2))));
    band_freqs = find(f1 >= band(1) & f1 <= band(2));
    pow = -trapz(f1(band_freqs), log10(pxx(band_freqs))) / (band(2) - band(1));
end

function [ freq ] = PeakFreq(f1, pxx, band)
    band_freqs_index = find(f1 >= band(1) & f1 <= band(2));
    [~, maxValIndex] = max(pxx(band_freqs_index));
    freqs = f1(band_freqs_index);
    freq = freqs(maxValIndex);
end
