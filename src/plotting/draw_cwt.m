function [] = draw_cwt(cfs,time,freq)
    args = {time, freq, cfs};
    surf(args{:},'edgecolor','none');
    view(0,90);
    axis xy;
    %axis tight;
    ax = gca;
    ax.XAxis.Visible = 'off'; % remove x-axis
    %ax.YTick = 0:10:200;
    
    shading interp; 
    %colormap(parula(256));
    colormap(jet);
    colorbar('off');
    %h = colorbar;
    %h.Label.String = 'z-score';
    ylabel('Frequency (Hz)');
end