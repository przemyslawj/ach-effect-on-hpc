function [] = showTraces( x, fs, figure_title, channelTable, shift )

if nargin < 3
    figure_title = '';
end

times = (1:size(x,2)) / fs;

vars = struct();
vars.shift = shift;
vars.x = x;
vars.times = times;
vars.fs = fs;
vars.lengthSeconds = size(x,2) / fs;
vars.max_timewindow = 4;

rec_len_sec = vars.max_timewindow * 0.5;
vars.channelTable = channelTable;
[~, vars.channelOrder] = sort(channelTable.channel_name);

figure('Name', figure_title,...
    'Position', [100 100 1200 800]);
nfigure = length(findobj('type','figure'));
plotTraces(0, rec_len_sec, vars);

sliderHandle = uicontrol('Style', 'slider', ...
          'Position', [10 20 800 20],...
          'SliderStep', [1 / vars.lengthSeconds / vars.max_timewindow,...
                         10 / vars.lengthSeconds / vars.max_timewindow],...
          'Tag', sprintf('timeSlider_%d', nfigure)); 
set(sliderHandle,'Callback',{@sliderCallback, vars});

zoomSliderHandle = uicontrol('Style', 'slider', ...
          'Position', [10 50 100 20],...
          'SliderStep', [0.1, 0.3],...
          'Tag', sprintf('zoomSlider_%d', nfigure));
zoomSliderHandle.Value = 0.5;
set(zoomSliderHandle,'Callback',{@sliderCallback, vars});
end

function [] = plotTraces(time_start_sec, rec_len_sec, vars)
    rec_end_i = min(size(vars.x,2), ceil((time_start_sec + rec_len_sec) * vars.fs));
    indecies = floor(time_start_sec * vars.fs + 1) : rec_end_i;
    nchan = size(vars.x, 1);
    minY = max(-10, min(vars.x));
    for i = 1:nchan
        ch = vars.channelOrder(i);
        trace = vars.x(ch,:);
        shifted_trace = (nchan - i) * vars.shift + trace;
        plot(vars.times(indecies), shifted_trace(indecies));
        hold on;
    end
    hold off;
    maxY = max(shifted_trace);
    %ylim([minY maxY]);
    grid on;
    grid minor;
    %legend(fliplr(vars.channelTable.channel_name));
    %legend(vars.channelTable.channel_name(vars.channelOrder));
    labels = strcat(...
        num2str(vars.channelTable.channel(vars.channelOrder)),...
        repmat('-', size(vars.channelOrder)),...
        vars.channelTable.channel_name(vars.channelOrder));
    legend(labels);
end

function [] = sliderCallback(hObject, evt, vars)
    updatePlotTraces(hObject, vars)
end

function [] = updatePlotTraces(hObject, vars)
    tag = get(hObject,'Tag');
    tag_parts = split(tag,'_');

    hZoom = findobj('Tag', ['zoomSlider_', tag_parts{2}]);
    rec_len_sec = hZoom.Value * vars.max_timewindow;
    hTime = findobj('Tag', ['timeSlider_', tag_parts{2}]);
    start_time_sec = (vars.lengthSeconds - rec_len_sec) * hTime.Value;
    plotTraces(start_time_sec, rec_len_sec, vars);
end
