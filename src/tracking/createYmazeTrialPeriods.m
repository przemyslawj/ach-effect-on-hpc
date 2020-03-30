function [trialPeriods] = createYmazeTrialPeriods(dataArray, time_mouse_arrived, fs)
    
key_position_names = { 'StartZone',  'FirstArm', 'Junction', 'SecondArm', ...
    'GoalZone', 'AfterReachedReward_10sec', 'AfterConsumedReward_10_sec', ...
    'BeforeGoalZone', 'DuringStim', 'DuringMaze', 'StimLast10sec', 'Total'};

key_positions_sample = round(time_mouse_arrived.sec * fs);
key_positions_sample(1) = max(key_positions_sample(1), 1);
if key_positions_sample(1) > size(dataArray, 2)
    return
end

consuming_reward_sec = 10;
after_reward_sample = (time_mouse_arrived.sec(end) + consuming_reward_sec) * fs;
interval_indecies = [...
    key_positions_sample(:), ... 
    [key_positions_sample(2:end); after_reward_sample]];  % AfterReachedReward_10sec

recording_end_sample = size(dataArray,2);
interval_indecies = [
    interval_indecies;
    [(after_reward_sample), (after_reward_sample + round(10 * fs))];  % AfterConsumerdReward_10_sec
    [key_positions_sample(1) key_positions_sample(end - 1)];  % BeforeGoalZone
    [key_positions_sample(end-1) recording_end_sample]; % DuringStim
    [key_positions_sample(1) key_positions_sample(end)];  % DuringMaze
    [max(recording_end_sample - 10 * fs, key_positions_sample(end-1)) recording_end_sample];  % StimLast10sec
    [key_positions_sample(1) recording_end_sample]  % Total
];

trialPeriods = table(interval_indecies(:,1),...
                     interval_indecies(:,2),...
                     key_position_names');
                
trialPeriods.Properties.VariableNames = {'starts', 'ends', 'stage_desc'};
end