function [] = sliderCallback(hObject, evt, vars)
    slider_pos = get(hObject,'Value');
    start_time_sec = (vars.lengthSeconds - 3) * slider_pos;
    plotData(vars.fs, vars.data, vars.filtered, vars.spikeFiltered, ...
             vars.ripple_starts, vars.ripple_ends, ...,
             vars.spikeTimes, vars.spikeThreshold, start_time_sec);
end

