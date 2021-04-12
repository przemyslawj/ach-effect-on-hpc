Code analysing data and generating figures for publication: *Impaired spatial learning and suppression of sharp wave ripples by cholinergic activation at the goal location*, https://doi.org/10.7554/eLife.65998.

The code is written in Matlab and R.

# Prerequisites
* install Matlab Waivelet Toolbox
* run src/compile.m

# Analysis of Y-maze learning
The learning curves were assessed with R notebook located in src/analysis/R/ymaze_learning.Rmd.

# Extracting Sharp-wave ripples (SWRs) and power spectrum density (PSD)
First, The sharp wave ripples and power spectrum density were extracted with Matlab using individual trial recordings and relevant periods of the trials. The results were saved to a .csv file and then analysed in R scripts.
For the first step data recorded with WHISPER system was analysed using src/analysis/matlab/state_psd_welch.m, data recorded with Igor, using src/analysis/matlab/state_psd_urethane_welch.m.

Further statistical analysis on the SWRs and PSD was assessed with R notebooks in src/analysis/R/ directory. Trials recorded during sleep: sleep_ripples.Rmd and sleep_ripples.Rmd; during ymaze task: ymaze_psd.Rmd and ymaze_ripples.Rmd; urethane recordings were analysed with: state_psd_urethane.Rmd.

# Plotting ephys recordings
## src/plotting/PlotElectrodeSignal.m
Shows raw signal of a single channel and the detected SWRs.
To run set in the script variables:
* *path* - set the top directory with data, for the given day
* *channelList* - vector of channels to display. Default is to show all channels.
  Alternatively use a csv file listing valid channel numbers for each animal
  using src/signal/findSelectedChannels.m function.

Plots from the top:
1) 80-250 Hz filtered signal with the detected ripples marked by black
horizontal bars
2) 80-250 Hz filtered signal squared - ripples are detected when the plot
crosses set threshold
3) Raw signal
Top scroll controls zoom, bottom scroll shifts the time of the plot.

## src/plotting/PlotMultipleElectrodeSignal.m
Shows raw signal on all of the electrodes.
To run set params as for PlotElectrodeSignal.m

## src/plotting/showSpectrogram.m
Shows the raw signal, SWR-filtered signal and spectrogram of the raw signal.
To run set params:
* *path*
* *channelList*
* *secondOffset* -- time of the recording when the plot is started
* *lengthSeconds* -- duration of the plot
