function [] = timewindowCallback(hObject, evt, vars)
rec_len_sec = get(hObject,'Value') * 4;
tag = get(hObject,'Tag');
tag_parts = split(tag,'_');

h = findobj('Tag', ['timeSlider_', tag_parts{2}]);
h.UserData.rec_len_sec = rec_len_sec;
start_time_sec = (vars.lengthSeconds - 2) * h.UserData.slider_pos;
plotData(vars.fs, vars.data, vars.filtered, ...
         vars.ripple_starts, vars.ripple_ends, vars.normalizedSquaredSignal,...
         start_time_sec,...
         rec_len_sec);
end

