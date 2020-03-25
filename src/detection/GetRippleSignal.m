function [ripple_detection_signal] = GetRippleSignal(filtered, frequency)
%windowLength = round(frequency/23);
% Window length adjusted on frequency 1025 Hz
windowLength = frequency/1025 * 11;
% make window of odd length
windowLength = floor(windowLength / 2) * 2 + 1; 

squaredSignal = filtered.^2;
window = ones(windowLength,1)/windowLength;
ripple_detection_signal = MovingAverageFilter(window, squaredSignal);
end