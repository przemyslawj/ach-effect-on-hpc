function [animals_dat] = summarize_ymaze_psd_dir(datarootdir, output_file)

if nargin < 2
    datarootdir = '/media/prez/DATA/Prez/y-maze/2018-02-20';
    output_file = 'psd_table_ymaze_2018-02-20_cwt.csv';
end
datarootdir
signalpath = [datarootdir filesep 'signal'];

binfiles = dir([signalpath filesep '*.bin']);

animals_dat = struct();
ntrials = 10;
nchans = 32;
psd_table = table();


for i = 1:numel(binfiles)
    filename_parts  = strsplit(binfiles(i).name, '_');
    animal = filename_parts{1};
    dstr = datestr(binfiles(1).date, 'yyyy-mm-dd');
    experiment = [dstr '_' filename_parts{1} '_' filename_parts{2} '_' filename_parts{3}];
        
    if strcmp(animal, 'Nigel')
        channels = [17];
    else
        channels = [15];
    end
    
    % Set up data for animal
    if ~isfield(animals_dat, animal)
        channel_std = zeros(ntrials, nchans);
        animal_dat = struct(...
            'ntrials', 0, ...
            'mpsd', zeros(37,11),...
            'channel_std', channel_std);
    else
        animal_dat = getfield(animals_dat, animal);
    end

    
    psd_out = ymaze_trial_psd(datarootdir, binfiles(i), channels);
    if isempty(psd_out.pow)
        continue
    end
    
    nrows = size(psd_out.pow, 1);
    psd_out.pow.experiment = cellstr(repmat(experiment, nrows, 1));
    psd_out.pow.trial = cellstr(repmat(filename_parts{3}, nrows, 1));
    psd_out.pow.animal = cellstr(repmat(animal, nrows, 1));
    psd_out.pow.date = cellstr(repmat(dstr, nrows, 1));
    animal_dat.mpsd = animal_dat.mpsd + psd_out.psd_xx;
    animal_dat.psd_all_bands = psd_out.psd_all_bands;
    animal_dat.ntrials = animal_dat.ntrials + 1;
    
    animals_dat = setfield(animals_dat, animal, animal_dat);
    psd_table = [psd_table; psd_out.pow];
end

writetable(psd_table, output_file);

end
