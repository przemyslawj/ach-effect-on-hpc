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

% filter invalid pos
tracking_dat_filtered = tracking_dat(tracking_dat.x >= 0,:);
movie_fs = 15;
tracking_dat_filtered.time_sec = tracking_dat_filtered.frame / movie_fs - secondOffset;

reward_percent = round(max(tracking_dat_filtered.total_percent) - 10);
key_positions_percent = [0 45 90 120 150 reward_percent];

key_positions_sec = zeros(size(key_positions_percent));
for i = 1:numel(key_positions_percent)
    j = find(tracking_dat_filtered.total_percent > key_positions_percent(i), 1);
    if isempty(j) || j == 0
        j = numel(tracking_dat_filtered.time_sec);
    end
    key_positions_sec(i) = tracking_dat_filtered.time_sec(j);
end
key_positions_sec(end) = min(key_positions_sec(end), max(tracking_dat_filtered.frame) / movie_fs);

time_mouse_arrived = table(key_positions_percent', key_positions_sec',...
    'VariableNames', {'percent', 'sec'}, ...
    'RowNames', {'Start', 'FirstArm', 'Junction', 'SecondArm', 'Goal', 'End'});

end

