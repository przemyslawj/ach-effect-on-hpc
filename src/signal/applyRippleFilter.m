function [filtered] = applyRippleFilter(data, channelTable, fs) 
    nyquist = fs / 2;
    filterOrder = 4;
    filterRipple = 20;
    
    filtered = zeros(size(data));
    emg_scaler = 20;
    laser_scaler = 5;
    if ~isempty(channelTable.location)
        for chan = 1:size(data,1)
            loc = channelTable.location{chan};
            x = data(chan,:);
            x = x - mean(x);
            if strcmp(loc, 'EMG') 
                filtered(chan,:) = x / max(x) / emg_scaler;
            elseif strcmp(loc, 'Laser')
                filtered(chan,:) = x / laser_scaler;
            else
                passband = [80 250];
                %passband = [100 250];
                %if strcmp(loc, 'CA3')
                %    passband = [150 250];
                %end
                [b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);
                filtered(chan,:) = filtfilt(b, a, x);
            end
        end
    end
end