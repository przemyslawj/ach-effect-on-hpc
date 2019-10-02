function [laserPeriods] = findLaserPeriods(laserSignal, max_indexdiff)
    laserSignal = abs([0 laserSignal 0]);
    laserOnThreshold = 100;
    laserOnIndex = find(laserSignal > laserOnThreshold);
    laserTurnedOff = laserOnIndex(diff([laserOnIndex numel(laserSignal) + max_indexdiff + 1]) > max_indexdiff);
    laserTurnedOn = laserOnIndex(diff([-max_indexdiff-1 laserOnIndex]) > max_indexdiff);
    
    laserOnPeriods = table();
    if numel(laserTurnedOn) > 0 
        laserOnPeriods = array2table([laserTurnedOn' laserTurnedOff'],...
            'VariableNames', {'start', 'end'});
        laserOnPeriods.laserOn = ones(numel(laserTurnedOn), 1);
    end
    
    laserOffPeriods = table();
    if numel(laserTurnedOn) > 0 
        laserOffPeriods = array2table(...
            [[0; laserTurnedOff'] ...
             [laserTurnedOn'; numel(laserSignal)]],...
            'VariableNames', {'start', 'end'});
        laserOffPeriods.laserOn = zeros(numel(laserTurnedOn) + 1, 1);
    end
    
    laserPeriods = [laserOffPeriods; laserOnPeriods];
    if size(laserPeriods, 1) > 0
        laserPeriods = sortrows(laserPeriods, 1);
    end
end