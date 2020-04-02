library(plyr)
library(dplyr)
library(purrr)
  
bin.ripples = function(ripple.df, bin_dur = 5) {
  df = ripple.df %>%
    dplyr::mutate(bin=ceiling(peak_t / bin_dur),
                  nripples=1)
  if (nrow(ripple.df) == 0) {
    return(ripple.df)
  }
  
  max_bin = max(df$bin)
  
  df.stage.template = df %>%
    dplyr::mutate(max.bin=max_bin) %>%
    ddply(.(animal, date, state, trial_id, channelLocation, channelName), 
          plyr::summarise,
          bin=1:max.bin[1])
  
  df.template = map_df(levels(df$stage_desc), ~ dplyr::mutate(df.stage.template, stage_desc = .x))
  df.template$stage_desc= as.factor(df.template$stage_desc)
  df.template$stage_desc = factor(df.template$stage_desc, levels=c('before_stim','stim','after_stim'))
  df.template = dplyr::mutate(df.template, laserOn=ifelse(stage_desc=='stim', 1, 0))
    df.template$laserOn = as.factor(df.template$laserOn)
  
  binned.ripple.df = df %>%
    right_join(df.template, by = c("animal", "date", "state", "channelName", "channelLocation", "trial_id", "stage_desc", "laserOn", "bin")) %>%
    replace_na(list(nripples=0)) %>% 
    dplyr::mutate(abs_bin = max_bin * (as.numeric(stage_desc) - 1) + bin) %>% 
    dplyr::group_by(trial_id, abs_bin, laserOn, channelLocation) %>%
    dplyr::summarise(nripples.total = sum(nripples, na.rm=TRUE),
              incidence = nripples.total / bin_dur) %>% 
    dplyr::mutate(bin_sec = abs_bin * bin_dur - bin_dur/2) 
  
  return(binned.ripple.df)
}