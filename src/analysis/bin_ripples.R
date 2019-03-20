bin.ripples = function(ripple.df, min_time, max_time, bin_dur = 5) {
  ripple.df = ripple.df %>%
    mutate(bin=ceiling(peak_t / bin_dur), bin_sec=bin * bin_dur - bin_dur/2)
  if (nrow(ripple.df) == 0) {
    return(ripple.df)
  }
  
  df = ripple.df %>% 
    mutate(nripples=1)
  
  max_bin = max(df$bin)
  min_bin = min(df$bin)
  nbins = max_bin - min_bin + 1
  trials = unique(df$trial_id)
  binned.template = data.frame(trial_id=rep(trials, each=nbins),
                               bin=rep(min_bin:max_bin, length(trials)),
                               nripples=rep(0, length(trials)))
  binned.ripple.df = df %>%
    right_join(binned.template, by=c('trial_id', 'bin')) %>%
    mutate(trial_id2=trial_id) %>%
    separate(trial_id2, c('cond', 'trial'), sep='_')
  binned.ripple.df$cond = as.factor(binned.ripple.df$cond)
  binned.ripple.df$trial = as.integer(binned.ripple.df$trial)
  
  binned.ripple.df = binned.ripple.df %>%
    group_by(trial_id, bin) %>%
    summarise(nripples.total = sum(nripples.x, na.rm=TRUE) + sum(nripples.y, na.rm=TRUE),
              incidence = nripples.total / bin_dur) %>%
    mutate(bin_sec = bin * bin_dur - bin_dur/2)
  
  return(binned.ripple.df)
}