library(pracma)
library(purrr)
library(reticulate)
library(stringr)

reticulate::use_condaenv('scipy')
np = reticulate::import('numpy', as='np')
fooof = reticulate::import('fooof')
fooof.analysis = reticulate::import('fooof.analysis')
fooof.synth = reticulate::import('fooof.synth')

calc.bands.pow = function(sample.psd.df, freq.bands.df, freq.range=c(3.5, 80), show.fit=FALSE, fit.knee=TRUE) {
  background.mode = 'fixed'
  if (fit.knee) {
    background.mode = 'knee'
  }
  fm = fooof$FOOOF(background_mode = background.mode, max_n_peaks = 8, verbose = TRUE,
                   peak_width_limits = c(1, 30))
  spectra.mat = np_array(c(sample.psd.df$band_power))
  freqs = np_array(sample.psd.df$band_start_freq)
  fm$fit(freqs = freqs, 
         power_spectrum = spectra.mat,
         freq_range = freq.range)
  res = list()
  for (band_i in 1:nrow(freq.bands.df)) {
    band.row = fooof.analysis$get_band_peak(fm$peak_params_, 
                                             band_def=c(freq.bands.df$xmin[band_i], freq.bands.df$xmax[band_i]), 
                                             ret_one = TRUE)
    #band.indecies = which(fm$freqs >= freq.bands.df$xmin[band_i] & fm$freqs <= freq.bands.df$xmax[band_i])
    #peak.freq = fm$freqs[band.indecies][which.max(fm$power_spectrum[band.indecies])]
    #band.rows = fooof.analysis$get_band_peak(fm$peak_params_, 
    #                                        band_def=c(freq.bands.df$xmin[band_i], freq.bands.df$xmax[band_i]), 
    #                                        ret_one = FALSE)
    #band.rows = matrix(band.rows, ncol=3)
    peak.freq = band.row[1]
    band_name = stringr::str_replace(freq.bands.df$band[band_i], ' ', '_')
    res[[paste('peak', band_name, sep='_')]] = peak.freq #weighted.mean(band.rows[,1], band.rows[,2])
    res[[paste('pow', band_name, sep='_')]] = band.row[2]
    res[[paste('width', band_name, sep='_')]] = band.row[3] #max(band.rows[,3])
  }
  res$background_offset = fm$background_params_[1]
  if (length(fm$background_params_) == 3) {
    res$background_knee = fm$background_params_[2]
    res$background_slope = fm$background_params_[3]
  } else {
    res$background_slope = fm$background_params_[2]
  }
  res$background_auc = trapz(log10(fm$freqs), fm[['_bg_fit']])
  res$fit_error = fm$error_
  #y = 10^offset * (1/(knee + x^exp))
  if (show.fit) {
    fm$print_results()
    matplot(fm$freqs, cbind(fm$power_spectrum, fm$fooofed_spectrum_, fm[['_bg_fit']]), 
            type='l',lty=1:5, log = "x")
    res$fm = fm
  }
  return(res)
}

test.psd.fitting = function() {
  #sample.trial.id = unique(psd.molten$trial_id)[1]
  sample.trial.id = '2019-11-29_GD_GD_sleep_2_g0_0'
  sample.psd.df = filter(psd.molten, 
                         trial_id == sample.trial.id & 
                           stage_desc == 'before_stim' & laserOn == 0 & 
                           #stage_desc == 'stim' & laserOn == 1 & 
                           channelLocation == 'CA1') %>%
    dplyr::arrange(band_start_freq)
  fres = calc.bands.pow(sample.psd.df, freq.bands.df, 
                        freq.range = c(3.5, 80),
                        show.fit = TRUE)
  fres
}

fit.background.df = function(psd.molten, freq.bands.df, freq.range = c(3.5, 80), fit.knee=TRUE) {
  psd.entries = psd.molten %>% 
    dplyr::select(animal, trial_id, stage_desc, laserOn, channelLocation) %>%
    dplyr::distinct()
  
  fit.psd.df = map_dfr(1:nrow(psd.entries), ~ {
    entry.df = psd.molten %>%
      filter(.data$trial_id == psd.entries$trial_id[.x] &
             .data$stage_desc == psd.entries$stage_desc[.x] & 
             .data$laserOn == psd.entries$laserOn[.x] & 
             .data$channelLocation == psd.entries$channelLocation[.x])
    #lowbands.pow.list = calc.bands.pow(entry.df, filter(freq.bands.df, band == 'theta'), 
    #                                   freq.range = c(1, 35), show.fit = F)
    #names(lowbands.pow.list) = paste('low', names(lowbands.pow.list),sep='_')
    highbands.pow.list = calc.bands.pow(entry.df, 
                                        freq.bands.df, #filter(freq.bands.df, band != 'theta'), 
                                        freq.range = freq.range, 
                                        show.fit = FALSE,
                                        fit.knee = fit.knee)
    return(cbind(psd.entries[.x,], #lowbands.pow.list, 
                 highbands.pow.list))
  })
  fit.psd.df
}

gen.mean.psd = function(param.group) {
  background.vals = fooof.synth$gen_background(
      fm$freqs, 
      c(mean(param.group$background_offset),
        mean(param.group$background_slope)))
  peak.theta = mean(param.group$peak_theta, na.rm=T)
  amp.theta = mean(param.group$pow_theta, na.rm=T)
  width.theta = mean(param.group$width_theta, na.rm=T)
  peak.sgamma = mean(param.group$peak_slow_gamma, na.rm=T)
  amp.sgamma = mean(param.group$pow_slow_gamma, na.rm=T)
  width.sgamma = mean(param.group$width_slow_gamma, na.rm=T)
  peak.mgamma = mean(param.group$peak_med_gamma, na.rm=T)
  amp.mgamma = mean(param.group$pow_med_gamma, na.rm=T)
  width.mgamma = mean(param.group$width_med_gamma, na.rm=T)
  peak.vals = fooof.synth$gen_peaks(fm$freqs, 
                                    gauss_params=np_array(c(peak.theta, amp.theta, width.theta,
                                                            peak.sgamma, amp.sgamma, width.sgamma,
                                                            peak.mgamma, amp.mgamma, width.mgamma)))
  yvals = peak.vals + background.vals
  return(yvals)
}

