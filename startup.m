thisDir = fileparts(mfilename('fullpath'));
beneDir = fullfile(thisDir, 'benedikt');

% Add current directory and its sub-directories to path
addpath(genpathGit(thisDir));

% Make sure that Benedikt's repository doesn't shadow anything
benePath = genpath(beneDir);

rmpath(benePath);
addpath(benePath, '-end');

SynEM.setup();

% Mark as ready
global PIPELINE_READY;
PIPELINE_READY = true;
clear;
