function [keep_sample] = excludeEMGNoisePeriods(emgsignal, crop_window_length, threshold)

if nargin < 3
    threshold = 2;
end

is_kept = ones(size(emgsignal));
sig_std = std(emgsignal);
diff_vals = abs(diff(emgsignal));

noise_peaks = find(diff_vals > sig_std * threshold);

halfwindow_len = int32(crop_window_length / 2);

for i = 1:numel(noise_peaks)
   peak_i = noise_peaks(i);
   noise_indecies = max(1, peak_i - halfwindow_len) : ...
                    min(numel(emgsignal), peak_i + halfwindow_len);
     
   is_kept(noise_indecies) = 0;
end

keep_sample = find(is_kept);
end