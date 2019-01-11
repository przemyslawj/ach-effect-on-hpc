function [] = showTraces( x, fs, figure_title, channelList )

if nargin < 3
    figure_title = '';
end

times = (1:size(x,2)) / fs;

shift = 5 * median(std(x(channelList,:),[],2));

vars = struct();
vars.shift = shift;
vars.x = x;
vars.times = times;
vars.fs = fs;
vars.lengthSeconds = size(x,2) / fs;

vars.rec_len_sec = 2;
vars.channelList = channelList;

figure('Name', figure_title);
plotTraces(0, vars);
sliderHandle = uicontrol('Style', 'slider', ...
          'Position', [10 20 500 20]); 
set(sliderHandle,'Callback',{@tracesSliderCallback, vars});
end

function [] = plotTraces(time_start_sec, vars)
    rec_end_i = min(size(vars.x,2), ceil((time_start_sec + vars.rec_len_sec) * vars.fs));
    indecies = floor(time_start_sec * vars.fs + 1) : rec_end_i;
    nchan = numel(vars.channelList);
    minY = max(-10, min(vars.x(vars.channelList,:)));
    i = 1;
    for ch = vars.channelList'
        trace = vars.x(ch,:);
        shifted_trace = (nchan - i) * vars.shift + trace;
        plot(vars.times(indecies), shifted_trace(indecies));
        hold on;
        i = i + 1;
    end
    hold off;
    maxY = max(shifted_trace);
    %ylim([minY maxY]);
    grid on;
    grid minor;
    legend(num2str(fliplr(vars.channelList)));
end

function [] = tracesSliderCallback(hObject, evt, vars) 
    slider_pos = get(hObject,'Value');
    start_time_sec = (vars.lengthSeconds - vars.rec_len_sec) * slider_pos;
    plotTraces(start_time_sec, vars)
end
