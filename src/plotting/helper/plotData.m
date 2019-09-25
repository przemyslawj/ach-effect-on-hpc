function [] = plotData(fs, raw, filtered, ...
                       ripple_starts, ripple_ends, normalizedSquaredSignal,...
                       start_time_sec, rec_len_sec)
    start_index = round(start_time_sec * fs) + 1;
    end_index = min(start_index + rec_len_sec * fs, length(filtered));
    indecies = start_index : end_index; 
    times = indecies / fs;
    
    subplot(3,1,1);
    plotSWR(times, filtered(indecies), fs, ripple_starts, ripple_ends);
    %ylim([-0.05 0.05]);

    subplot(3,1,2);
    plot(times, normalizedSquaredSignal(indecies));
    %ylim([-0.5 6]);
    
    subplot(3,1,3);
    plot(times, raw(indecies));
    %ylim([-0.3 0.5]);
    
end

