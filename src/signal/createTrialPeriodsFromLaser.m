function [trialPeriodsTable] = createTrialPeriodsFromLaser(laserPeriods, durationPeriod)

% Check laser off and on alternate and finish with laser off
nperiods = size(laserPeriods, 1);
expectedLaser = zeros(nperiods, 1);
expectedLaser(2:2:nperiods) = 1;
assert(all(expectedLaser == laserPeriods.laserOn))

nlaserOn = sum(laserPeriods.laserOn==1);

before_laser_indecies = 1:2:(nperiods-1);
laser_indecies = 2:2:nperiods;
after_laser_indecies = 3:2:nperiods;


starts = [
    max(laserPeriods.start(before_laser_indecies),...
        laserPeriods.end(before_laser_indecies) - durationPeriod);
    laserPeriods.start(laser_indecies);
    laserPeriods.start(after_laser_indecies) ];

ends = [
    laserPeriods.end(before_laser_indecies);
    min(laserPeriods.start(laser_indecies) + durationPeriod,...
        laserPeriods.end(laser_indecies));
    min(laserPeriods.start(after_laser_indecies) + durationPeriod,...
        laserPeriods.end(after_laser_indecies)) ];

stage_desc = [
    repmat({'before_stim'}, nlaserOn, 1);
    repmat({'stim'}, nlaserOn, 1);
    repmat({'after_stim'}, nlaserOn, 1) ];

laserOn = [
    zeros(nlaserOn, 1);
    ones(nlaserOn, 1);
    zeros(nlaserOn, 1)];


trialPeriodsTable = table(starts, ends, stage_desc, laserOn);
trialPeriodsTable = sortrows(trialPeriodsTable, 'starts');
end
