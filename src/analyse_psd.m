ntrials = animal_dat.ntrials;
psd_pow = animal_dat.pow_theta;
channels_std = animals_dat.Albert.channel_std;
channels_std = channels_std(mean(channels_std,2) > 0,:);
key_position_names = animal_dat.position_names;
fs = size(psd_pow, 3);
nchan = size(psd_pow, 2);

all_names = cell(1, 2 * numel(key_position_names));
for i = 1:numel(key_position_names)
    name = key_position_names{i};
    all_names{2 * i - 1} = strcat(name, '_mean');
    all_names{2 * i} = strcat(name, '_zscore');
end

pow_theta_summary = zeros(nchan, numel(all_names));
for channel = 1:nchan
    channel_std = channels_std(:,channel);
    keep_trials = find(channel_std < 0.2);
    
    channel_pow = squeeze(psd_pow(keep_trials, channel, :));
    means_indecies = 1:2:numel(all_names);
    means = mean(channel_pow, 'omitnan');
    pow_theta_summary(channel,means_indecies) = means;
    std_indecies = 2:2:numel(all_names);
    pow_theta_summary(channel,std_indecies) = std(channel_pow, 'omitnan') ./ means;
end

pow_theta_summary = array2table(pow_theta_summary);
pow_theta_summary.Properties.VariableNames = all_names;
row_names = cellfun(@(x) num2str(x), num2cell(0:31),'UniformOutput',false);
pow_theta_summary.Properties.RowNames = row_names;