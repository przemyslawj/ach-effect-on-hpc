%% setup results
datarootdir = '/mnt/DATA/Clara/baseline/';

sleeping_trials_fpath = [datarootdir filesep 'trials.csv'];
expstable = readtable(sleeping_trials_fpath, 'ReadVariableNames', true);
%expstable = expstable(strcmp(expstable.animal, animal_code),:);

nchans = 4;

nexp = size(expstable, 1);
result_table = table();

all_ripples = table();

result_table.animal = repmat('XX', nexp * 3 * nchans, 1);
result_table.state = repmat('xxxxx', nexp * 3 * nchans, 1);
result_table.date = repmat('2018-01-01', nexp * 3 * nchans, 1);
result_table.trial = zeros(nexp * 3 * nchans, 1);
result_table.laserOn = zeros(nexp * 3 * nchans, 1);
result_table.channel = zeros(nexp * 3 * nchans, 1);
all_bands = exp(0.7:0.05:5.3);
result_table.all_psd_xx = zeros(nexp * 3 * nchans, numel(all_bands));

result_table.pow_slow = zeros(nexp * 3 * nchans, 1);
result_table.pow_theta = zeros(nexp * 3 * nchans, 1);
result_table.peak_theta = zeros(nexp * 3 * nchans, 1);
result_table.pow_slow_gamma = zeros(nexp * 3  * nchans, 1);
result_table.peak_slow_gamma = zeros(nexp * 3  * nchans, 1);
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
    animal_code = expstable.animal{i};
    state = expstable.state{i};
    electrodes_file = '/mnt/DATA/Clara/ymaze/valid_electrodes.csv';
    electrodes = readtable(electrodes_file);
    channelList = electrodes(strcmp(electrodes.animal, animal_code),:).channel;
    nchans = numel(channelList);
    
    filename_prefix = [animal_code '_trial_' num2str(expstable.trial(i)) '_'];
    dateddir = datestr(expstable.date(i), 'yyyy-mm-dd');
    binfile = dir([datarootdir filesep dateddir filesep 'signal' filesep filename_prefix '*.bin']);
    fprintf('Processing file: %s\n', binfile.name);
    meta = ReadMeta(binfile.name, binfile.folder);

    break_sec = 0;
    duration_sec = 20;
    total_sec = str2double(meta.fileTimeSecs);
    rec_times = [expstable.startLaserSec(i) - duration_sec, expstable.startLaserSec(i);...
                 expstable.startLaserSec(i) + break_sec, expstable.startLaserSec(i) + break_sec + duration_sec;...
                 expstable.stopLaserSec(i) + break_sec , min(total_sec - 2, expstable.stopLaserSec(i) + break_sec + duration_sec)];

    
    %rec_times = [expstable.startLaserSec(i) - duration_sec - break_sec, expstable.startLaserSec(i) - break_sec;...
    %         expstable.startLaserSec(i) + break_sec, expstable.startLaserSec(i) + break_sec + duration_sec;...
    %         expstable.stopLaserSec(i) + break_sec, min(total_sec - 2, expstable.stopLaserSec(i) + break_sec + duration_sec)];
    %rec_times = [10, 20;...
    %             30, 40;
    %             expstable.stopLaserSec(i) + break_sec, min(total_sec - 2, expstable.stopLaserSec(i) + break_sec + duration_sec)];
    


    dataArray = ReadSGLXData(meta, 0, rec_times(end));
    dataArray = dataArray(channelList, :);

    fs = 1250;
    time = (1:size(dataArray,2)) / fs;
    dataArrayOrg = dataArray;
    dataArray = downsample(dataArray', round(meta.nSamp / fs))';
    nchans = size(dataArray, 1);

    %dataArray = filter50Hz(dataArray, fs);
    channel_std = std(dataArray, [], 2);
    downfs = 625;
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

    slow = [3 7];
    theta = [7 12];
    above_theta = [12 25];
    slow_gamma = [25 45];
    med_gamma = [60 80];
    fast_gamma = [80 150];

    for channel = 1:nchans
        
        [ripples, sd, normalizedSquaredSignal] = MyFindRipples(time', filtered(channel,:)', ...
                     'frequency', fs, ...
                     'thresholds', [2 5 0.01],...
                     'durations', [10 30 350]);
        
        for rec_times_i = 1:size(rec_times, 1)
            sec_start = rec_times(rec_times_i, 1);
            sec_end = rec_times(rec_times_i, 2);
            sec_length = sec_end - sec_start;

            if sec_length <= 0
                entry_i = entry_i + 1;
                continue
            end

            result_table.channel(entry_i,:) = channelList(channel);
            result_table.date(entry_i,:) = dateddir;
            result_table.animal(entry_i,:) = animal_code;
            result_table.state(entry_i,:) = state;
            result_table.trial(entry_i) = expstable.trial(i);
            result_table.laserOn(entry_i) = mod(rec_times_i - 1, 2);
            
            if channel_std(channel) > 0.1
                entry_i = entry_i + 1;
                continue
            end

            psd_xx = downsampledDataArray(channel,sec_start * downfs + 1:sec_end * downfs);

            %[pxx, freqs] = pwelch(psd_xx, ...
            %    floor(downfs / 4), floor(downfs / 8), floor(downfs / 2), downfs);
            [cws, freqs] = cwt(psd_xx, 'amor', downfs);        
            pxx = median(abs(cws .^ 2), 2);
            freqs = fliplr(freqs')';
            pxx = fliplr(pxx')';

            result_table.dom_freq(entry_i) = PeakFreq(freqs, pxx, [slow(1), 45]);

            result_table.pow_slow(entry_i) = TotalBandPower(freqs, pxx, slow);
            result_table.pow_theta(entry_i) = TotalBandPower(freqs, pxx, theta);
            result_table.peak_theta(entry_i) = PeakFreq(freqs, pxx, theta);
            result_table.pow_slow_gamma(entry_i) = TotalBandPower(freqs, pxx, slow_gamma);
            result_table.peak_slow_gamma(entry_i) = PeakFreq(freqs, pxx, slow_gamma);
            result_table.pow_med_gamma(entry_i) = TotalBandPower(freqs, pxx, med_gamma);
            result_table.pow_fast_gamma(entry_i) = TotalBandPower(freqs, pxx, fast_gamma);
            result_table.pow_above_theta(entry_i) = TotalBandPower(freqs, pxx, above_theta);
            for j=1:(numel(all_bands)-1)
                result_table.all_psd_xx(entry_i, j) = TotalBandPower(freqs, pxx, [all_bands(j) all_bands(j+1)]);
            end

             
             section_ripples = ripples(ripples(:,2) >= sec_start & ripples(:,2) <= sec_end,:);
             nripples = size(section_ripples, 1);
             if ~isempty(section_ripples)
                result_table.has_ripples(entry_i) = 1;
                result_table.swr_incidence(entry_i) = size(section_ripples,1) / sec_length;
                result_table.ripple_peakpow(entry_i) = mean(section_ripples(:,4));
                result_table.ripple_dur(entry_i) = mean(section_ripples(:,3) - section_ripples(:,1));
                result_table.ripple_freq(entry_i) = mean(section_ripples(:,5));
               
                ripples_table = array2table(section_ripples);
                ripples_table.Properties.VariableNames = {'start_sec', 'peak_t', 'end_time', 'peakpow', 'peak_freq'};
                ripples_table.date = repmat(dateddir, nripples, 1);
                ripples_table.animal = repmat(animal_code, nripples, 1);
                ripples_table.trial = repmat(expstable.trial(i), nripples, 1);
                ripples_table.channel = repmat(channelList(channel), nripples, 1);
                ripples_table.laserOn = repmat(result_table.laserOn(entry_i), nripples, 1);
                ripples_table.state = repmat(state, nripples, 1);
 
                all_ripples = [all_ripples; ripples_table];
             end
             result_table.swr_powerhz(entry_i) = mean(abs(normalizedSquaredSignal));
            
             entry_i = entry_i + 1;
        end
       
    end
end

writetable(result_table, [datarootdir filesep 'psd_table_th4.csv']);
writetable(all_ripples, [datarootdir filesep 'ripples_th4.csv']);
            
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
