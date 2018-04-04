animals_dat3 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-09', 'psd_table_ymaze_2018-02-09_cwt.csv');
animals_dat4 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-10', 'psd_table_ymaze_2018-02-10_cwt.csv');
animals_dat5 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-11', 'psd_table_ymaze_2018-02-11_cwt.csv');
animals_dat1 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-07', 'psd_table_ymaze_2018-02-07_cwt.csv');
animals_dat2 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-08', 'psd_table_ymaze_2018-02-08_cwt.csv');
animals_dat6 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-12', 'psd_table_ymaze_2018-02-12_cwt.csv');
animals_dat7 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-13', 'psd_table_ymaze_2018-02-13_cwt.csv');
animals_dat8 = summarize_ymaze_psd_dir('/media/prez/DATA/Prez/y-maze/2018-02-20', 'psd_table_ymaze_2018-02-20_cwt.csv');

total_trials = animals_dat1.Albert.ntrials + animals_dat2.Albert.ntrials + ...
    animals_dat3.Albert.ntrials + animals_dat4.Albert.ntrials + animals_dat5.Albert.ntrials + ...
    animals_dat6.Albert.ntrials + animals_dat7.Albert.ntrials + animals_dat8.Albert.ntrials;
mpsd = animals_dat1.Albert.mpsd + animals_dat2.Albert.mpsd + animals_dat3.Albert.mpsd + animals_dat4.Albert.mpsd + ...
    animals_dat5.Albert.mpsd + animals_dat6.Albert.mpsd + animals_dat7.Albert.mpsd + animals_dat8.Albert.mpsd;
mpsd = mpsd / total_trials;
mpsdtable = array2table(mpsd);
mpsdtable.Properties.VariableNames = key_position_names;
mpsdtable.freqs_start = animals_dat1.Albert.psd_all_bands(1:end-1)';
mpsdtable.freqs_end = animals_dat1.Albert.psd_all_bands(2:end)';

writetable(mpsdtable, 'albert_mpsd.csv');


total_trials = animals_dat1.Nigel.ntrials + animals_dat2.Nigel.ntrials + animals_dat3.Nigel.ntrials + animals_dat4.Nigel.ntrials + ...
    animals_dat5.Nigel.ntrials + animals_dat6.Nigel.ntrials + animals_dat7.Nigel.ntrials + animals_dat8.Nigel.ntrials;
mpsd = animals_dat1.Nigel.mpsd + animals_dat2.Nigel.mpsd + animals_dat3.Nigel.mpsd + animals_dat4.Nigel.mpsd + ...
    animals_dat5.Nigel.mpsd + animals_dat6.Nigel.mpsd + animals_dat7.Nigel.mpsd + animals_dat8.Nigel.mpsd;
mpsd = mpsd / total_trials;
mpsdtable = array2table(mpsd);
mpsdtable.Properties.VariableNames = key_position_names;
mpsdtable.freqs_start = animals_dat1.Nigel.psd_all_bands(1:end-1)';
mpsdtable.freqs_end = animals_dat1.Nigel.psd_all_bands(2:end)';

writetable(mpsdtable, 'nigel_mpsd.csv');
