
bands = c(exp(seq(0.7, 5.30, by=0.05)))

melt.freq.bands = function(df, extra.id.vars=c('learning_day', 'dirname')) {
  default.id.vars = c('animal', 'date', 'trial_id', 'channelLocation','laserOn', 'stage_desc')
  
  id.vars = c(default.id.vars, extra.id.vars)
  psd.molten = df %>%
    select(id.vars, starts_with('all_psd_xx')) %>% 
    reshape2::melt(id.vars=id.vars, variable.name='band_name', value.name='band_power')  %>%
    dplyr::mutate(band_start_freq = bands[as.numeric(stringr::str_replace(band_name, 'all_psd_xx_', ''))])
}

calc.freq.band.ratio = function(psd.molten) {
  reshape2::dcast(psd.molten,
                  animal + date + trial_id + channelLocation + band_name + band_start_freq ~ stage_desc,
                  value.var = 'band_power') %>%
    dplyr::mutate(power.ratio.change=stim / before_stim)
}