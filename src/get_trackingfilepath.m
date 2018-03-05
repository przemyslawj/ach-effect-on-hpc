function [filepath] = get_trackingfilepath(datarootdir, binName)
trackingdir = [datarootdir filesep 'movie' filesep 'tracking'];
filename_parts = strsplit(binName, '_');
tracking_filename = ['processed_' strjoin(filename_parts(1:3), '_') '_positions.csv'];
filepath = [trackingdir filesep tracking_filename];
end

