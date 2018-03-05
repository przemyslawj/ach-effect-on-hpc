datarootdir = '/home/przemek/neurodata/LFP_data/y-maze/2018-02-13';
signalpath = [datarootdir filesep 'signal'];

binfiles = dir([signalpath filesep '*.bin']);

animals_dat = struct();
ntrials = 10;
nchans = 32;

for i = 1:numel(binfiles)
    filename_parts  = strsplit(binfiles(i).name, '_');
    animal = filename_parts{1};
    
    % Set up data for animal
    if ~isfield(animals_dat, animal)
        animal_pow_theta = zeros(ntrials, nchans, 10);
        channel_std = zeros(ntrials, nchans);
        animal_dat = struct(...
            'ntrials', 0, ...
            'channel_std', channel_std);
    else
        animal_dat = getfield(animals_dat, animal);
    end
        
    psd_out = ymaze_trial_psd(datarootdir, binfiles(i));
    pow_theta = psd_out.pow_theta;
    
    if isempty(pow_theta)
        continue
    end
    animal_dat.ntrials = animal_dat.ntrials + 1;
    animal_dat.position_names = pow_theta.Properties.VariableNames;
    animal_dat.pow_slow(animal_dat.ntrials, :, :) = table2array(psd_out.pow_slow);
    animal_dat.pow_theta(animal_dat.ntrials, :, :) = table2array(pow_theta);
    animal_dat.pow_above_theta(animal_dat.ntrials, :, :) = table2array(psd_out.pow_above_theta);
    animal_dat.slow_gamma(animal_dat.ntrials, :, :) = table2array(psd_out.pow_slow_gamma);
    animal_dat.channel_std(animal_dat.ntrials,:) = psd_out.channel_std;
    animals_dat = setfield(animals_dat, animal, animal_dat);
end

save('animals_dat_2018-02-13.mat', 'animals_dat');

