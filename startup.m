% Add current directory and its sub-directories to path
addpath(genpathGit(fileparts(mfilename('fullpath'))));

% Mark as ready
global PIPELINE_READY;
PIPELINE_READY = true;

