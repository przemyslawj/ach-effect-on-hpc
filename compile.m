% this function recompiles all .c scripts to create the apropriate .mex* files on your local machine
%NOTE: current directory must be the buzcode base directory

addpath(genpath('externalPackages'))
addpath(genpath('src'))

try
cd('externalPackages/FilterM/')
catch
    display('Please navigate to your main directory of swr detect')
end

% mex -O FilterX.c % compiles FilterM (faster filtering than filtfilt)
mex('CFLAGS="\$CFLAGS -std=c99"', 'FilterX.c') % the above line fails with newer compilers but this works
cd('../..')
