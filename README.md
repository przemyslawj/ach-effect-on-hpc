# Requirements for Matlab
* install Matlab Waivelet Toolbox

# Instructions
* run src/compile.m to use the code

# Plotting
## src/plotting/PlotElectrodeSignal.m
Shows raw signal of a single channel and the detected SWRs.
To run set in the script variables:
* *path* - set the top directory with data, for the given day
* *channelList* - vector of channels to display. Default is to show all channels.
  Alternatively use a csv file listing valid channel numbers for each animal
  using src/signal/findSelectedChannels.m function.

Plots from the top:
1) 100-250 Hz filtered signal with the detected ripples marked by black
horizontal bars
2) 100-250 Hz filtered signal squared - ripples are detected when the plot
crosses set threshold
3) Raw signal
4) 300 Hz High pass filtered signal for spike detection

Top scroll controls zoom, bottom scroll shifts the time of the plot.

## src/plotting/PlotMultipleElectrodeSignal.m
Shows raw signal on all of the electrodes.
To run set params as for PlotElectrodeSignal.m

## src/plotting/showSpectrogram.m
Shows raw signal, SWR-filtered signal and spectrogram of the raw signal.
To run set params:
* *path*
* *channelList*
* *secondOffset* -- time of the recording when the plot is started
* *lengthSeconds* -- duration of the plot
