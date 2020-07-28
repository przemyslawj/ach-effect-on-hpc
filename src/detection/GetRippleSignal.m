% Returns signal after applying moving average filter on the squared signal.
% Moving average windows is equal to duration of 10.4 ms
function [ripple_detection_signal] = GetRippleSignal(filtered, frequency)
%windowLength = round(frequency/23);
% Window length of 10.4 ms
windowLength = frequency/1025 * 13;
% make window of odd length
windowLength = floor(windowLength / 2) * 2 + 1; 

squaredSignal = filtered.^2;
window = ones(windowLength,1)/windowLength;
ripple_detection_signal = MovingAverageFilter(window, squaredSignal);
end