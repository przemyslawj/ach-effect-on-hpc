% Basic SWR detection based on a threshold, use MyFindRipples.m instead.

function [ timestamps, peak_starts, peak_endings ] = detectSWR(x, fs)

threshold = 3 * std(x);
%min_peak_dur = round(0.005 * fs);
min_peak_dur = 0;

timestamps = double(x > threshold);

timestamp_str = ['0' sprintf('%d', timestamps) '0'];
peak_starts = strfind(timestamp_str, '01');
peak_endings = strfind(timestamp_str, '10');
for i = 1:length(peak_starts)
    peak_dur = peak_endings(i) - peak_starts(i) + 1;
    if (peak_dur <= min_peak_dur)
        timestamps((peak_starts(i) - 1) : (peak_endings(i) - 1)) = ...
            zeros(1, peak_dur);
    end
end

end
