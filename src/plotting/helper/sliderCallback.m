function [] = sliderCallback(hObject, evt, vars)
    hObject.UserData.slider_pos = get(hObject,'Value');
	
    start_time_sec = (vars.lengthSeconds - 2) * hObject.UserData.slider_pos;
    plotData(vars.fs, vars.data, vars.filtered, ...
             vars.ripple_starts, vars.ripple_ends, vars.normalizedSquaredSignal,...
             start_time_sec,...
             hObject.UserData.rec_len_sec);
end

