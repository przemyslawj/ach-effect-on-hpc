%% setup results
datarootdir = '/media/prez/DATA/Prez/N&A_rest';
statefilepath = [datarootdir filesep 'state.csv'];
statetable = readtable(statefilepath, 'ReadVariableNames', true);
statetable = statetable(strcmp(statetable.animal, 'Albert'), :);
animalstate = 'rest';
statetable = statetable(strcmp(statetable.state, animalstate), :);
experiment_names = unique(statetable.experiment);

nexp = size(statetable, 1);
nnames = numel(experiment_names);
pxxlen = 70;
psd_on = zeros(nnames, pxxlen);
psd_off = zeros(nnames, pxxlen);
psd_off_adjacent = zeros(nnames, pxxlen);


%% process single experiments
for i = 1:nexp
    
    if strcmp(statetable.laser{i}, 'OFF') && statetable.time_sec(i) > 0
        continue;
    end
    
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
    
    [pxx, freqs] = ReadDataAndCalculatePSD(meta, psd_sec_start, psd_sec_end - psd_sec_start);
    
    exp_index = find(strcmp(experiment_names, statetable.experiment{i}));
    if strcmp(statetable.laser{i}, 'ON')
        psd_on(exp_index,1:pxxlen) = pxx(1:pxxlen);
    else
        psd_off(exp_index,1:pxxlen) = pxx(1:pxxlen);
        [pxx, freqs] = ReadDataAndCalculatePSD(meta, psd_sec_start - trial_duration, psd_sec_end - psd_sec_start);
        psd_off_adjacent(exp_index,1:pxxlen) = pxx(1:pxxlen);
    end

end

%% averages
freqs=freqs(1:pxxlen);
figure('Name', ['State: ' animalstate]);

subplot(2,2,1);
i=1;
PlotPSDExamples(freqs, psd_on(i,1:pxxlen), psd_off(i,1:pxxlen), psd_off_adjacent(i,1:pxxlen));
ylabel('LFP Power (db/Hz)')

subplot(2,2,2);
i=3;
PlotPSDExamples(freqs, psd_on(i,1:pxxlen), psd_off(i,1:pxxlen), psd_off_adjacent(i,1:pxxlen));

subplot(2,2,3);
%psd_diff = psd_on(:,1:pxxlen) ./ psd_off(:,1:pxxlen);
psd_diff = (log10(psd_on(:,1:pxxlen)) - log10(psd_off(:,1:pxxlen))) ./ log10(psd_on(:,1:pxxlen));
%psd_diff_adjacent = psd_off_adjacent(:,1:pxxlen) ./ psd_off(:,1:pxxlen);
psd_diff_adjacent = (log10(psd_off_adjacent(:,1:pxxlen)) - log10(psd_off(:,1:pxxlen))) ./ log10(psd_off(:,1:pxxlen));


plot(freqs, mean(psd_diff, 'omitnan'), 'b');
hold on;
for i = 1:size(psd_diff,1)
    lp = semilogx(freqs, psd_diff(i,:), 'black');
    lp.Color(4) = 0.3;
end
sd = std(psd_diff, 'omitnan');
%plot(freqs, mean(psd_diff, 'omitnan') + sd, 'r');
%plot(freqs, mean(psd_diff, 'omitnan') - sd, 'r');
xlim([0 50]);
xlabel('Frequency (Hz)');
ylabel('LFP power change Stimulation/Control (%)');
ylim([-0.3  0.3]);
xlim([0 50]);

hold off;
subplot(2,2,4);
plot(freqs, mean(psd_diff_adjacent, 'omitnan'), 'b');
hold on;
for i = 1:size(psd_diff,1)
    lp = semilogx(freqs, psd_diff_adjacent(i,:), 'black');
    lp.Color(4) = 0.3;
end
sd = std(psd_diff_adjacent, 'omitnan');
%plot(freqs, mean(psd_diff_adjacent, 'omitnan') + sd, 'r');
%plot(freqs, mean(psd_diff_adjacent, 'omitnan') - sd, 'r');
xlim([0 50]);
ylim([-0.2 0.2]);
xlabel('Frequency (Hz)')
ylabel('LFP power change Control1/Control2 (%)');

function [] = PlotPSDExamples(freqs, psd_on, psd_off, psd_off_adjacent)
    semilogx(freqs, log10(psd_off), 'Color','b');
    hold on;
    semilogx(freqs, log10(psd_off_adjacent), 'Color','b', 'LineStyle','--');
    semilogx(freqs, log10(psd_on), 'Color','r');
    xlim([0 50]);
    %ylim([-5 -3.2]);
end

function [pxx, freqs] = ReadDataAndCalculatePSD(meta, psd_sec_start, sec_length)
    dataArray = ReadSGLXData(meta, psd_sec_start, sec_length);

    fs=256;
    fs=1250;
    downsampledDataArray = downsample(dataArray', round(meta.nSamp / fs))';
    downsampledDataArray = filter50Hz(downsampledDataArray, fs);

    %% calculate PSD
    channel = 16;

    psd_x = downsampledDataArray(channel,:);
    if isempty(psd_x)
        return
    end

    [cws, freqs] = cwt(psd_x, 'morse', fs);
    pxx = mean(abs(cws .^ 2), 2);
    %[pxx, freqs] = pwelch(psd_x, ...
    %        floor(fs / 4), floor(fs / 4) - floor(fs/8), floor(fs / 8), fs);
    %[pxx, freqs] = pwelch(psd_x,floor(fs / 4), floor(fs / 4) - floor(fs/8), 1:100, fs);
    
end  