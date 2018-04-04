function [] = plotData(fs, raw, filtered, spikeFiltered, ...
                       ripple_starts, ripple_ends, normalizedSquaredSignal,...
                       spikeTimes, spike_threshold, start_time_sec)
    rec_len_sec = 2;
    start_index = round(start_time_sec * fs) + 1;
    end_index = min(start_index + rec_len_sec * fs, length(filtered));
    indecies = start_index : end_index; 
    times = indecies / fs;
    
    subplot(4,1,1);
    plotSWR(times, filtered(indecies), fs, ripple_starts, ripple_ends);
    ylim([-0.05 0.05]);

    subplot(4,1,2);
    plot(times, normalizedSquaredSignal(indecies));
    ylim([-0.5 6]);
    
    subplot(4,1,3);
    plot(times, raw(indecies));
    ylim([-0.3 0.5]);
    
    subplot(4,1,4);
    plot(times, spikeFiltered(indecies));
    ylim([-0.08 0.08]);
    
    % Plot spikes
    hold on;
    maxY = 0.07;
    for i = 1:length(spikeTimes)
        if spikeTimes(i) >= times(1) && spikeTimes(i) <= times(end)
            line([spikeTimes(i) - 0.005, spikeTimes(i) + 0.005], [maxY, maxY], 'Color', 'black', 'LineWidth', 6);
        end;
    end
    
    line([times(1), times(end)], [spike_threshold spike_threshold], 'Color', 'red');
    line([times(1), times(end)], -[spike_threshold spike_threshold], 'Color', 'red');
    hold off;
end

