function [ trialPeriods ] = extractTrialPeriodsFromLaser(...
    downsampledDataArray, laserChannelIdx, downfs, laserOnThreshold)
    lasertrial_duration_sec = 20;
    laserPeriods = [];
    if laserChannelIdx > 0
        laserSignal = downsampledDataArray(laserChannelIdx, :);
        laserPeriods = findLaserPeriods(laserSignal, downfs * 0.2, laserOnThreshold);
    end
    
    trialPeriods = [];
    if ~isempty(laserPeriods)
        trialPeriods = createTrialPeriodsFromLaser(laserPeriods, ...
            downfs * lasertrial_duration_sec);
    end
end