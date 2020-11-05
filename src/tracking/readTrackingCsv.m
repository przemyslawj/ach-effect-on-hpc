function [time_mouse_arrived] = readTrackingCsv(filepath, secondOffset)

if nargin < 2
    secondOffset = 0;
end

time_mouse_arrived = table();
if ~exist(filepath, 'file')
    warning('Tracking file %s not found\n', filepath)
    return
end
    
tracking_dat = readtable(filepath, 'ReadVariableNames', true, ...
    'PreserveVariableNames', true);

% filter invalid pos and positions with jumps
tracking_dat_filtered = tracking_dat(tracking_dat.x > 0 & tracking_dat.dist < 50, :);
tracking_dat_filtered.time_sec = tracking_dat_filtered.timestamp / 1000;

goal_percent = 150;
reward_percent = round(max(tracking_dat_filtered.total_percent) - 10);
key_positions_percent = [0 45 90 120 goal_percent reward_percent];

key_positions_sec = zeros(size(key_positions_percent));
key_positions_index = zeros(size(key_positions_percent));
for i = 1:numel(key_positions_percent)
    j = find(tracking_dat_filtered.total_percent > key_positions_percent(i), 1);
    if isempty(j) || j == 0
        j = numel(tracking_dat_filtered.time_sec);
    end
    key_positions_sec(i) = tracking_dat_filtered.time_sec(j);
    key_positions_index(i) = j;
end

arrived_end_index = key_positions_index(end);
goal_arm = tracking_dat_filtered.arm(arrived_end_index);
left_goal_index = find(...
        [ tracking_dat_filtered.total_percent(arrived_end_index+1:end); 
          goal_percent - 1 ]...
        < goal_percent ...
      & [ tracking_dat_filtered.arm(arrived_end_index+1:end); goal_arm ]...
        == goal_arm, 1) + arrived_end_index - 1;
key_positions_sec = [key_positions_sec tracking_dat_filtered.time_sec(left_goal_index)];
key_positions_sec = max(0, key_positions_sec - secondOffset);
key_positions_percent = [key_positions_percent goal_percent];
time_mouse_arrived = table(key_positions_percent', key_positions_sec',...
    'VariableNames', {'percent', 'sec'}, ...
    'RowNames', {'Start', 'FirstArm', 'Junction', 'SecondArm', 'Goal', 'End', 'LeftGoal'});

end

