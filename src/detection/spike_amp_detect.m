function [ spike_times, spike_indecies, thr] = spike_amp_detect( x, xf_detect )
%SPIKE_AMP_DETECT Finds spikes based on amplitude threshold for the recording.
% Threshold based on median estimation, spikes detected when positive or
% negative threshold is crossed.
%
% Noisy spikes are detected if they cross threshold set for the signal 
% bandpassed on lower frequencies.
%   Params
%   x - vector with the recorded signal
%   x_detect - vector with signal bandpassed for detection

% refractory period
vars = struct();
vars.fs = 20000;
vars.ref = round(0.005 * vars.fs);
vars.stdmin = 5;
vars.stdmax = 50;

% Bandpass x between 300 and 3000 for noise removal
filterOrder = 4;
filterRipple = 40;
nyquist = vars.fs / 2;
[b_sort, a_sort] = cheby2(filterOrder,filterRipple,[300 3000]/nyquist);
xf_sort = filtfilt(b_sort, a_sort, x);

noise_std_detect = median(abs(xf_detect))/0.6745;
% threshold for spike detection
thr = vars.stdmin * noise_std_detect;

noise_std_sort = median(abs(xf_sort))/0.6745;
% maximum threshold used for artifact removal.
thrmax = vars.stdmax * noise_std_sort;     

spike_times_i = zeros(1, ceil(length(x) / vars.ref));
nspikes = 0;
threshold_i = find(-xf_detect(1:end) > thr);
previous_spike_i = 0;
for i = 1:length(threshold_i)
    if threshold_i(i) >= previous_spike_i + vars.ref
        xaux_end = min(threshold_i(i) + floor(vars.ref/2) - 1, length(xf_detect));
        [~, maxamp_i] = max((xf_detect(threshold_i(i):xaux_end)));
        if maxamp_i == 1
            continue
        end
        nspikes = nspikes + 1;
        spike_times_i(nspikes) = maxamp_i + threshold_i(i) -1;
        previous_spike_i = spike_times_i(nspikes);
    end
end

%Eliminate artifacts
spikes_keep = zeros(size(spike_times_i));
for i = 1:nspikes
    if max(abs( xf_sort (spike_times_i(i)) )) <= thrmax
        spikes_keep(i) = 1;
    end
end

spike_times_i = spike_times_i(spikes_keep > 0);
%nnoisy = nspikes - length(spike_times_i);
%fprintf('Detected: %d spikes after removing %d noisy\n', length(spike_times_i), nnoisy);

spike_times = spike_times_i / vars.fs;
spike_indecies = spike_times_i;
end
