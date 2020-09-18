function [trialPeriods] = createYmazeTrialPeriods(dataArray, time_mouse_arrived, fs)
    
key_position_names = { 'StartZone',  'FirstArm', 'Junction', 'SecondArm', ...
    'GoalZone', 'AfterReachedReward_10sec', 'AfterConsumedReward_10_sec', ...
    'BeforeGoalZone', 'DuringStim', 'MazeStem', 'StimLast10sec', 'Total', ...
    'GoalZone20sec'};

key_positions_sample = round(time_mouse_arrived.sec * fs);
key_positions_sample(1) = max(key_positions_sample(1), 1);
if key_positions_sample(1) > size(dataArray, 2)
    return
end

consuming_reward_sec = 10;
after_reward_sample = (time_mouse_arrived.sec(end-1) + consuming_reward_sec) * fs;
interval_indecies = [...
    key_positions_sample(1:end-1), ... 
    [key_positions_sample(2:end-1); after_reward_sample]];  % AfterReachedReward_10sec

goal_zone_sample = key_positions_sample(end - 2);
ommit_last_sec = 2;
recording_end_sample = size(dataArray,2) - ommit_last_sec * fs;
left_goal_sample = time_mouse_arrived.sec(end) * fs;
interval_indecies = [
    interval_indecies;
    [(after_reward_sample), (after_reward_sample + round(10 * fs))];  % AfterConsumerdReward_10_sec
    [key_positions_sample(1) goal_zone_sample];  % BeforeGoalZone
    [key_positions_sample(end-2) left_goal_sample]; % DuringStim
    [key_positions_sample(2) goal_zone_sample];  % MazeStem
    [max(left_goal_sample - (10 - ommit_last_sec) * fs, goal_zone_sample) left_goal_sample];  % StimLast10sec
    [key_positions_sample(1) left_goal_sample];  % Total
    [goal_zone_sample (goal_zone_sample + 20 * fs)] % First 20 sec in goal zone
];

trialPeriods = table(interval_indecies(:,1),...
                     interval_indecies(:,2),...
                     key_position_names');
                
trialPeriods.Properties.VariableNames = {'starts', 'ends', 'stage_desc'};
end