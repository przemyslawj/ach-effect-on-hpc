function [filepath] = get_trackingfilepath(dateddir, binName)
trackingdir = [dateddir filesep 'movie' filesep 'tracking'];
dirparts = strsplit(dateddir, filesep);
dirparts(cellfun('isempty',dirparts)) = [];
datestr = dirparts{numel(dirparts)};
filename_parts = strsplit(binName, '_');
tracking_filename = ['processed_' strjoin([datestr filename_parts(1:3)], '_') '_positions.csv'];
filepath = [trackingdir filesep tracking_filename];
end

