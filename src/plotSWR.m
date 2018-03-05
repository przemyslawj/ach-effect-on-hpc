function [] = plotSWR(x, y, fs, swr_starts, swr_endings)
    plot(x, y);
    hold on;
    
    if isempty(swr_starts)
        swrThreshold = 3 * std(y);
        line([x(1) max(x)], [swrThreshold swrThreshold], 'Color', 'r');
        [~, swr_starts_i, swr_endings_i] = detectSWR(y,fs);
        swr_starts = x(swr_starts_i);
        swr_endings = x(swr_endings_i);
    end
    maxY = min(0.045, max(y) * 1.4);
    for i = 1:length(swr_starts)
        if swr_starts(i) >= x(1) && swr_endings(i) <= x(end)
            line([swr_starts(i), swr_endings(i)], [maxY, maxY], 'Color', 'black', 'LineWidth', 6);
        end
    end
    hold off;
end

