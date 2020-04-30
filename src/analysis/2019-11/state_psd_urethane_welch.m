%% setup results
rootdir = '/mnt/DATA/chat_ripples/urethane/';
experiment = 'scopolamine';
datarootdir = [ rootdir experiment filesep ];

txtfiles = dir([datarootdir filesep '*.txt']);

nchans = 1;
nexp = size(txtfiles, 1);
result_table = table();

all_ripples = table();
all_bands = exp(0.1:0.035:5.6) - 1;

result_table.animal = repmat('XXXXXX', nexp * 3 * nchans, 1);

result_table.row = zeros(nexp * 3 * nchans, 1);
result_table.trial = zeros(nexp * 3 * nchans, 1);
result_table.laserOn = zeros(nexp * 3 * nchans, 1);
result_table.channel = zeros(nexp * 3 * nchans, 1);
result_table.all_psd_xx = zeros(nexp * 3 * nchans, numel(all_bands));

result_table.pow_slow = zeros(nexp * 3 * nchans, 1);
result_table.pow_theta = zeros(nexp * 3 * nchans, 1);
result_table.peak_theta = zeros(nexp * 3 * nchans, 1);

result_table.pow_slow_gamma = zeros(nexp * 3  * nchans, 1);
result_table.peak_slow_gamma = zeros(nexp * 3 * nchans, 1);
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

    fname = [datarootdir filesep txtfiles(i).name]
    animal = txtfiles(i).name(1:6);
    datatable = readtable(fname);
    dataArray = datatable{:,1};
    dataArray = dataArray';
    lengthSeconds = size(dataArray, 2) / fs;

    rec_times = [0 15;
                 15 30;
                 45 60; ];    
%    rec_times = [0 60;
%                 60 120;
%                 120 180; ];

    time = (1:size(dataArray,2)) / fs;

    channel_std = std(dataArray, [], 2);

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
    
    stage_descs = {'before_stim', 'stim', 'after_stim'};

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
            result_table.stage_desc(entry_i) = stage_descs(rec_times_i);

            psd_xx = dataArray(channel,...
                sec_start * fs + 1:min(size(dataArray,2), sec_end * fs));

            [pxx, freqs] = pwelch(psd_xx, ...
                hanning(fs/2), floor(fs / 4), all_bands, fs);

            [maxVal, maxValIndex] = max(pxx);
            result_table.dom_freq(entry_i) = freqs(maxValIndex);

            result_table.pow_slow(entry_i) = CalcBandPower(fs, pxx, slow);
            result_table.pow_theta(entry_i) = CalcBandPower(fs, pxx, theta);
            result_table.peak_theta(entry_i) = PeakFreq(freqs, pxx, theta);
            result_table.pow_slow_gamma(entry_i) = CalcBandPower(fs, pxx, slow_gamma);
            result_table.peak_slow_gamma(entry_i) = PeakFreq(freqs, pxx, slow_gamma);
            result_table.pow_med_gamma(entry_i) = CalcBandPower(fs, pxx, med_gamma);
            result_table.pow_fast_gamma(entry_i) = CalcBandPower(fs, pxx, fast_gamma);
            result_table.pow_above_theta(entry_i) = CalcBandPower(fs, pxx, above_theta);
            for j=1:(numel(all_bands)-1)
                result_table.all_psd_xx(entry_i, j) = pxx(j);
            end

             filtered_start_i = sec_start * fs + 1;
             filtered_end_i = min(sec_end * fs + 1, size(filtered,2));
             [section_ripples, sd, normalizedSquaredSignal] = MyFindRipples(...
                     (time(filtered_start_i : filtered_end_i))', ...
                     (filtered(channel, filtered_start_i : filtered_end_i))', ...
                     'frequency', fs, ...
                     'thresholds', [2 4.0 0.01], ...
                     'durations', [10 20 300], ....
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
                ripples_table.stage_desc = repmat(stage_descs(rec_times_i), nripples, 1);
 
                all_ripples = [all_ripples; ripples_table];
             end
             result_table.swr_powerhz(entry_i) = mean(abs(normalizedSquaredSignal));
            
             entry_i = entry_i + 1;
        end
       
    end
end

outdir = [rootdir filesep 'trial_results'];
writetable(result_table, [outdir filesep 'welch_psd_table_urethane_' experiment '.csv']);
writetable(all_ripples, [outdir filesep 'ripples_urethane_' experiment '.csv']);
            

function [ pow ] = CalcBandPower(fs, signal, band)
    pow = bandpower(signal, fs, band);
end

function [ freq ] = PeakFreq(f1, pxx, band)
    band_freqs_index = find(f1 >= band(1) & f1 <= band(2));
    [~, maxValIndex] = max(pxx(band_freqs_index));
    freqs = f1(band_freqs_index);
    freq = freqs(maxValIndex);
end
