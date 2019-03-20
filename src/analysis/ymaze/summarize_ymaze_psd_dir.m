function [animals_dat] = summarize_ymaze_psd_dir(datarootdir, output_file)

if nargin < 2
    datarootdir = '/mnt/DATA/Clara/ymaze/2018-08-17';
    output_file = 'psd_table_ymaze_2018-08-17_cwt.csv';
end
datarootdir
signalpath = [datarootdir filesep 'signal'];

binfiles = dir([signalpath filesep '*.bin']);

animals_dat = struct();
ntrials = 10;
nchans = 32;
psd_table = table();
rewarded_arms = readtable('/mnt/DATA/Clara/ymaze/rewarded_arm.csv');
electrodes_file = '/mnt/DATA/Clara/ymaze/selected_electrodes.csv';
invalid_trials = readtable('/mnt/DATA/Clara/ymaze/invalid_trials.csv');
electrodes = readtable(electrodes_file);

for i = 1:numel(binfiles)
    filename_parts  = strsplit(binfiles(i).name, '_');
    animal = filename_parts{1};
    dstr = datestr(binfiles(1).date, 'yyyy-mm-dd');
    trial_num = str2double(filename_parts{3});
    experiment = [dstr '_' filename_parts{1} '_' filename_parts{2} '_' filename_parts{3}];
    
    trial_id_str = [dstr '_' animal '_' trial_num];
    if sum(strcmp(invalid_trials.trial_id, trial_id_str)) > 0
        sprintf('Invalid trial skipped: %s', trial_id_str)
        continue
    end

    
    channels = electrodes(strcmp(electrodes.animal, animal),:).channel;
    nchans = numel(channels);
    
    % Set up data for animal
    if ~isfield(animals_dat, animal)
        channel_std = zeros(ntrials, nchans);
        animal_dat = struct(...
            'ntrials', 0, ...
            'all_psd', [],...
            'channel_std', channel_std,...
            'ripples', []);
    else
        animal_dat = getfield(animals_dat, animal);
    end

    
    psd_out = ymaze_trial_psd(datarootdir, binfiles(i), channels);
    if isempty(psd_out.pow) || any(any(psd_out.all_psd_xx > 0)) == 0
        continue
    end
    
    rewarded_arm = rewarded_arms{strcmp(rewarded_arms.animal, animal),2};
    tracking_filepath = get_trackingfilepath(datarootdir, binfiles(i).name);
    tracking_dat = readtable(tracking_filepath, 'ReadVariableNames', true);
    tracking_dat_filtered = tracking_dat(tracking_dat.x >= 0,:);
    success_trial = tracking_dat_filtered.arm(end) == rewarded_arm && ...
        tracking_dat_filtered.total_percent(end) >= 150;
    
    nrows = size(psd_out.pow, 1);
    psd_out.pow.experiment = cellstr(repmat(experiment, nrows, 1));
    psd_out.pow.trial = cellstr(repmat(filename_parts{3}, nrows, 1));
    psd_out.pow.animal = cellstr(repmat(animal, nrows, 1));
    psd_out.pow.date = cellstr(repmat(dstr, nrows, 1));
    psd_out.pow.success = repmat(success_trial, nrows, 1);
    if (success_trial)
        animal_dat.all_psd = cat(3, animal_dat.all_psd, psd_out.all_psd_xx);
        animal_dat.psd_all_bands = psd_out.psd_all_bands;
        animal_dat.ntrials = animal_dat.ntrials + 1;
        nripples = size(psd_out.ripples, 1);
        if nripples > 0
            psd_out.ripples.trial = repmat(trial_num, nripples, 1);
            psd_out.ripples.animal = repmat(animal, nripples, 1);
            psd_out.ripples.date = repmat(dstr, nripples, 1);
            animal_dat.ripples = cat(1, animal_dat.ripples, psd_out.ripples);
        end
    end
    %animal_dat.ntrials = animal_dat.ntrials + psd_out.ntrials;
    
    animals_dat = setfield(animals_dat, animal, animal_dat);
    psd_table = [psd_table; psd_out.pow];
end

writetable(psd_table, output_file);

end

