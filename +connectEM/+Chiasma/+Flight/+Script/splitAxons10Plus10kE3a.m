% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.Util.loadFlightPaths;
import connectEM.Chiasma.Util.loadTaskIds;

clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a');

% task specifics
taskGenId = '20171108T105021';

taskGenFile = fullfile( ...
    chiasmaDir, 'taskGeneration', ...
    sprintf('%s_taskGeneration.mat', taskGenId));
taskIdFile = fullfile( ...
    chiasmaDir, 'taskAnswers', ...
    sprintf('%s_flightTaskIDs.txt', taskGenId));
nmlDir = fullfile( ...
    chiasmaDir, 'taskAnswers', ...
    sprintf('%s_flightTasks', taskGenId));

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

taskGenData = load(taskGenFile);

% task IDs
taskIds = loadTaskIds(taskIdFile);
taskIds = taskIds.id;

ff = loadFlightPaths(param, nmlDir);

% sanity checks
assert(size(taskGenData.exits, 1) == numel(taskIds));