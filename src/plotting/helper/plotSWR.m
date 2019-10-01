function [] = plotSWR(x, y_norm, fs, swr_starts, swr_endings)
    plot(x, y_norm);
    hold on;
    %maxY = min(0.045, max(y) * 1.4);
    maxY = max(y_norm) * 1.10;
    for i = 1:length(swr_starts)
        if swr_starts(i) >= x(1) && swr_endings(i) <= x(end)
            line([swr_starts(i), swr_endings(i)], [maxY, maxY], 'Color', 'black', 'LineWidth', 6);
        end
    end
    hold off;
end

