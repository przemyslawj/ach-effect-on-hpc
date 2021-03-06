
#bands = c(exp(seq(0.1, 5.6, by=0.07)))
bands = c(exp(seq(0.1, 5.6, by=0.035)))

melt.freq.bands = function(df, 
                           extra.id.vars=c('learning_day', 'dirname'), 
                           psd.vars.prefix='all_psd_xx_',
                           freq.bands = bands) {
  default.id.vars = c('animal', 'date', 'trial_id', 'channelLocation','laserOn', 'stage_desc')
  
  id.vars = c(default.id.vars, extra.id.vars)
  psd.molten = df %>%
    select(id.vars, starts_with(psd.vars.prefix)) %>% 
    reshape2::melt(id.vars=id.vars, variable.name='band_name', value.name='band_power')  %>% 
    dplyr::mutate(band_start_freq = freq.bands[as.numeric(stringr::str_replace(band_name, psd.vars.prefix, ''))])
}

calc.freq.band.ratio = function(psd.molten) {
  reshape2::dcast(psd.molten,
                  animal + date + trial_id + channelLocation + band_name + band_start_freq + state + exp ~ stage_desc,
                  value.var = 'band_power') %>%
    dplyr::mutate(power.ratio.change=stim / before_stim)
}

prepare.ymaze.df = function(df) {
  df$animal = as.factor(df$animal)
  df$channel = as.factor(df$channel)
  df$laserOn = as.factor(df$laserOn)
    
  df = df %>%
    dplyr::mutate(trial_id=paste(date, animal, trial ,sep='_'))
  
  df$trial_id = as.factor(df$trial_id)
  df = df %>%
    dplyr::group_by(date, animal) %>%
    dplyr::mutate(trial_no = map_dbl(strsplit(file_name, '_'), ~ as.numeric(.x[3])) ) 
  
  return(df)
}
