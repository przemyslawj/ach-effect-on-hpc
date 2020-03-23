path = '/mnt/DATA/Clara/urethane/light/';
[binName, path] = uigetfile('*.txt', 'LFP file', path);
fprintf('Processing file: %s\n', binName);

subplots = 0;

secondOffset = 0;

org_fs = 1000;
nChans = 1;
channelList = 1;

datatable = readtable([path filesep binName]);
dataArray = datatable{:,1};
dataArray = dataArray';

fs = 500;
dataArray = downsample(dataArray', round(org_fs / fs))';


lengthSeconds = 0.5;
startSec = 30.4;
timeIndecies = (startSec * fs) : ((startSec + lengthSeconds) * fs);
%timeIndecies = 1:size(dataArray, 2);
%timeIndecies = 1:fs*(140);

%dataArray = filter50Hz(dataArray, fs);

%% Filter LFP for ripples
filtered = zeros(size(dataArray));
passband = [80 200];
nyquist = fs / 2;
filterOrder = 4;
filterRipple = 20;
[b, a] = cheby2(filterOrder,filterRipple,passband/nyquist);

for chan = 1:nChans
    %filtered(i,:) = FiltFiltM(b, a, dataArray(i,:));
    % Use MATLAB built-in function
    filtered(chan,:) = filtfilt(b, a, dataArray(chan,:));
end

for i = 1:numel(channelList)
    channel = channelList(i);
    figName = [binName '-channel-' num2str(channel)];
    % Raw signal
    figure('Name', [figName ' 1']);
    subplot(4,1,1);

    timepoints = (1:size(dataArray,2)) / fs;
    plot(timepoints(timeIndecies), dataArray(channel,timeIndecies));
    %ylim([-0.5 0.5]);
    %ylim([-1000 1000]);
    xlabel('Time (sec)');  
    
    % Ripple signal
    subplot(4,1,2);
    %plot(timepoints, filtered(channel,:));
    [ripples,sd, normalizedSquaredSignal] = MyFindRipples(timepoints', filtered', ...
                                 'frequency', fs, ...
                                 'thresholds', [2 4 0.01],...
                                 'durations', [10 30 350]);
    ripple_starts = [];
    ripple_ends = [];
    if ~isempty(ripples)                               
        ripple_starts = ripples(:,1);
        ripple_ends = ripples(:,3);
    end
    plotSWR(timepoints(timeIndecies), filtered(timeIndecies), fs, ripple_starts, ripple_ends);
    ylim([-12 20]);
    % Spectrogram signal
    x = dataArray(channel,:);

   
    if subplots == 0
        figure('Name', [figName ' 2']);
    end
    %[wt, wfreqs]=cwt(x, 'morse', fs, 'ExtendSignal', true, ...
    %    'VoicePerOctave', 4, 'WaveletParameters', [3 120]);
    [wt, wfreqs]=cwt(x, 'amor', fs, 'ExtendSignal', true);
    wt_pow = abs(wt).^2;
    wt_zscored = z_score(wt_pow);
    high_freqs = find(wfreqs < 200);    
    subplot(4,1,3);
    draw_cwt(wt_zscored(high_freqs,timeIndecies), timepoints(timeIndecies), wfreqs(high_freqs));
    
    %[wt, wfreqs]=cwt(x, 'morse', fs, 'ExtendSignal', true, ...
    %                 'VoicePerOctave', 4, 'WaveletParameters', [3 40]);
    [wt, wfreqs]=cwt(x, 'amor', fs, 'ExtendSignal', true);
    wt_pow = abs(wt).^2;
    low_freqs = find(wfreqs <= 12);
    subplot(4,1,4);
    wt_zscored = z_score(wt_pow);
    draw_cwt(wt_zscored(low_freqs,timeIndecies), timepoints(timeIndecies), wfreqs(low_freqs));
    
    
%     figure('Name', ['Pwelch' '-channel-' num2str(channel)]);
%     [pxx, freqs] = pwelch(dataArray(channel,1:fs*15), ...
%                 floor(fs / 4), floor(fs / 8), floor(fs / 2), fs);
%     plot(freqs, 10*log10(pxx))
%     
%     figure('Name', ['CWT' '-channel-' num2str(channel)]);
%     [wt, wfreqs]=cwt(dataArray(channel,fs*1:fs*15), 'amor', fs);
%     wt_pow = abs(wt).^2;
%     plot(wfreqs, median(wt_pow, 2))
end

function A = z_score(cfs)
    A = cfs - mean(cfs, 2);
    A = bsxfun(@rdivide, A, std(A, [], 2));
    maxZscore = ones(size(A)) * 6;
    minZscore = ones(size(A)) * -2;
    A = bsxfun(@min, A, maxZscore);
    A = bsxfun(@max, A, minZscore);
end

function [] = draw_cwt(cfs,time,freq)
    args = {time, freq, cfs};
    surf(args{:},'edgecolor','none');
    view(0,90);
    axis xy;
    axis tight;
    ax = gca;
    %ax.XAxis.Visible = 'off'; % remove x-axis
    %ax.YTick = 0:10:200;
    
    shading interp; 
    %colorbar('off');
    h = colorbar;
    h.Label.String = 'z-score';
    ylabel('Frequency (Hz)');
    h.Limits= [-2 6];
    caxis([-2 6]);
    colormap(jet);
end


