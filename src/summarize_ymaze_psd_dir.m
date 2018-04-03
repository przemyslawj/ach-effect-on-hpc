function [] = summarize_ymaze_psd_dir(datarootdir, output_file)

if nargin < 2
    datarootdir = '/media/prez/DATA/Prez/y-maze/2018-02-20';
    output_file = 'psd_table_ymaze_2018-02-20_cwt.csv';
end

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
    
    % Set up data for animal
    if ~isfield(animals_dat, animal)
        channel_std = zeros(ntrials, nchans);
        animal_dat = struct(...
            'ntrials', 0, ...
            'channel_std', channel_std);
    else
        animal_dat = getfield(animals_dat, animal);
    end
    
    if strcmp(animal, 'Nigel')
        channels = [17];
    else
        channels = [15];
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
    psd_table = [psd_table; psd_out.pow];
end

writetable(psd_table, output_file);

end
