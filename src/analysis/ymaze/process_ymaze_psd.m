key_position_names = { 'StartZone',  'FirstArm', 'Junction', 'SecondArm', ...
    'GoalZone', 'AfterReachedReward_10sec', 'AfterConsumedReward_10_sec', ...
    'BeforeGoalZone', 'DuringStim', 'DuringMaze', 'HomeCageLast10sec', 'Total'};

animals_dat = [...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-15', 'psd_table_ymaze_2018-08-15_cwt.csv')...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-16', 'psd_table_ymaze_2018-08-16_cwt.csv')...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-17', 'psd_table_ymaze_2018-08-17_cwt.csv')...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-18', 'psd_table_ymaze_2018-08-18_cwt.csv')...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-19', 'psd_table_ymaze_2018-08-19_cwt.csv')...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-20', 'psd_table_ymaze_2018-08-20_cwt.csv')...
    summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-21', 'psd_table_ymaze_2018-08-21_cwt.csv')...
    %summarize_ymaze_psd_dir('/mnt/DATA/Clara/ymaze/2018-08-22', 'psd_table_ymaze_2018-08-22_cwt.csv')...
 ];

%%
animal_names = fields(animals_dat);
for field_i = 1:numel(animal_names)
    animal_name = animal_names{field_i};
    all_psd = [];
    total_trials = 0;
    for i = 1:size(animals_dat,2)
        animal_dat = animals_dat(i).(animal_name);
        if ~isempty(animal_dat.all_psd)
            all_psd = cat(3, all_psd, animal_dat.all_psd);
            total_trials = total_trials + animal_dat.ntrials;
        end
    end
    %wide_psd = wide_psd ./ total_trials;
    %wide_psdtable = psd2table(wide_psd, animal_dat.psd_wide_bands, key_position_names);
    
    mean_psd = median(all_psd, 3, 'omitnan');
    all_psdtable = psd2table(mean_psd, animal_dat.psd_all_bands, key_position_names);
    writetable(all_psdtable, ['all_psd_' animal_name '.csv']);
    
    std_psd = std(all_psd, [], 3, 'omitnan');
    std_psdtable = psd2table(std_psd, animal_dat.psd_all_bands, key_position_names);
    writetable(std_psdtable, ['std_psd_' animal_name '.csv']);
end

%%


function psd_table = psd2table(psd_dat, freq_bands, key_position_names)
    psd_table = array2table(psd_dat);
    psd_table.Properties.VariableNames = key_position_names;
    psd_table.freqs_start = freq_bands(1:end-1)';
    psd_table.freqs_end = freq_bands(2:end)';
end
