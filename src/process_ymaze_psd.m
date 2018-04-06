animals_dat1 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-07', 'psd_table_ymaze_2018-02-07_cwt.csv');
animals_dat2 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-08', 'psd_table_ymaze_2018-02-08_cwt.csv');
animals_dat3 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-09', 'psd_table_ymaze_2018-02-09_cwt.csv');
animals_dat4 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-10', 'psd_table_ymaze_2018-02-10_cwt.csv');
animals_dat5 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-11', 'psd_table_ymaze_2018-02-11_cwt.csv');
animals_dat6 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-12', 'psd_table_ymaze_2018-02-12_cwt.csv');
animals_dat7 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-13', 'psd_table_ymaze_2018-02-13_cwt.csv');
animals_dat8 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-20', 'psd_table_ymaze_2018-02-20_cwt.csv');

total_trials = animals_dat1.Albert.ntrials + animals_dat2.Albert.ntrials + ...
    animals_dat3.Albert.ntrials + animals_dat4.Albert.ntrials + animals_dat5.Albert.ntrials + ...
    animals_dat6.Albert.ntrials + animals_dat7.Albert.ntrials + animals_dat8.Albert.ntrials;
wide_psd = animals_dat1.Albert.wide_psd + animals_dat2.Albert.wide_psd + animals_dat3.Albert.wide_psd + animals_dat4.Albert.wide_psd + ...
    animals_dat5.Albert.wide_psd + animals_dat6.Albert.wide_psd + animals_dat7.Albert.wide_psd + animals_dat8.Albert.wide_psd;
wide_psd = wide_psd / total_trials;
wide_psdtable = psd2table(wide_psd, animals_dat1.Albert.psd_wide_bands, key_position_names);
all_psd = animals_dat1.Albert.all_psd + animals_dat2.Albert.all_psd + animals_dat3.Albert.all_psd + animals_dat4.Albert.all_psd + ...
    animals_dat5.Albert.all_psd + animals_dat6.Albert.all_psd + animals_dat7.Albert.all_psd + animals_dat8.Albert.all_psd;
all_psd = all_psd / total_trials;
all_psdtable = psd2table(all_psd, animals_dat1.Albert.psd_all_bands, key_position_names);

writetable(wide_psdtable, 'wide_psd_albert.csv');
writetable(all_psdtable, 'all_psd_albert.csv');


total_trials = animals_dat1.Nigel.ntrials + animals_dat2.Nigel.ntrials + animals_dat3.Nigel.ntrials + animals_dat4.Nigel.ntrials + ...
    animals_dat5.Nigel.ntrials + animals_dat6.Nigel.ntrials + animals_dat7.Nigel.ntrials + animals_dat8.Nigel.ntrials;
wide_psd = animals_dat1.Nigel.wide_psd + animals_dat2.Nigel.wide_psd + animals_dat3.Nigel.wide_psd + animals_dat4.Nigel.wide_psd + ...
    animals_dat5.Nigel.wide_psd + animals_dat6.Nigel.wide_psd + animals_dat7.Nigel.wide_psd + animals_dat8.Nigel.wide_psd;
wide_psd = wide_psd / total_trials;
wide_psdtable = psd2table(wide_psd, animals_dat1.Nigel.psd_wide_bands, key_position_names);
all_psd = animals_dat1.Nigel.all_psd + animals_dat2.Nigel.all_psd + animals_dat3.Nigel.all_psd + animals_dat4.Nigel.all_psd + ...
    animals_dat5.Nigel.all_psd + animals_dat6.Nigel.all_psd + animals_dat7.Nigel.all_psd + animals_dat8.Nigel.all_psd;
all_psd = all_psd / total_trials;
all_psdtable = psd2table(all_psd, animals_dat1.Nigel.psd_all_bands, key_position_names);

writetable(wide_psdtable, 'nigel_wide_psd.csv');
writetable(all_psdtable, 'nigel_all_psd.csv');

function psd_table = psd2table(psd_dat, freq_bands, key_position_names)
    psd_table = array2table(psd_dat);
    psd_table.Properties.VariableNames = key_position_names;
    psd_table.freqs_start = freq_bands(1:end-1)';
    psd_table.freqs_end = freq_bands(2:end)';
end