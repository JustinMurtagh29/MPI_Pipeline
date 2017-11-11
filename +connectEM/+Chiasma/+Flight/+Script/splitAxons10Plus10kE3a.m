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

splitFile = fullfile( ...
    chiasmaDir, sprintf('%s_splitData.mat', taskGenId));
splitFileBuild = false;

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% build or load split file
if splitFileBuild || ~exist(splitFile, 'file')
    % task definitions
    taskGenData = load(taskGenFile);
    taskDefs = taskGenData.taskDefs;
    exits = taskGenData.exits;

    % chiasmata
    chiasmata = load(taskGenData.info.param.chiasmataFile);
    axonFile = chiasmata.info.param.axonFile;
    chiasmata = chiasmata.chiasmata;

    % task IDs
    taskIds = loadTaskIds(taskIdFile);
    taskIds = taskIds.id;

    % flight paths
    flights = loadFlightPaths(param, nmlDir);

    % build split file
    splitData = connectEM.Chiasma.Flight.prepareSplit( ...
        chiasmata, taskDefs, exits, taskIds,flights);
    splitData.axonFile = axonFile;
    splitData.info = info;
    
    % save
    Util.saveStruct(splitFile, splitData);
else
    splitData = load(splitFile);
end