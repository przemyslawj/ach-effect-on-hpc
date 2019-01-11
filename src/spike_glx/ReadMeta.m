function [meta] = ReadMeta(binName, path)

    % Create the matching metafile name
    [dumPath,name,dumExt] = fileparts(binName);
    metaName = strcat(name, '.meta');

    % Parse ini file into cell entries C{1}{i} = C{2}{i}
    fid = fopen(fullfile(path, metaName), 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
    C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
    fclose(fid);

    % New empty struct
    meta = struct();

    % Convert each cell entry into a struct entry
    for i = 1:length(C{1})
        tag = C{1}{i};
        if tag(1) == '~'
            % remake tag excluding first character
            tag = sprintf('%s', tag(2:end));
        end
        meta = setfield(meta, tag, C{2}{i});
    end
    
    %meta.typeThis = 'spikegl';
    %meta.nSamp = str2double(meta.niSampRate) * 0.8;
    %meta.fileTimeSecs = num2str(str2double(meta.fileTimeSecs) * 5/4);
    % Assumes correcting sampling set
    meta.nSamp = str2double(meta.niSampRate);
    meta.fileTimeSecs = num2str(meta.fileTimeSecs);
    meta.nChans = str2double(meta.nSavedChans);
    %meta.nChans = 32;
    meta.binName = binName;
    meta.path = path;
end % ReadMeta




