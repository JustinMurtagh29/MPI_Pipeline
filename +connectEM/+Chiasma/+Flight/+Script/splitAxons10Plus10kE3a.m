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

outFile = fullfile( ...
    chiasmaDir, sprintf('%s_splitAxons.mat', taskGenId));

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% build or load split file
if splitFileBuild || ~exist(splitFile, 'file')
    splitData = buildSplitData(taskGenFile, taskIdFile, nmlDir);
    splitData.info = info;
    
    Util.saveStruct(splitFile, splitData);
    clear splitData;
end

%% running chiasma splitting
axonFile = load(splitFile, 'axonFile');
axonFile = axonFile.axonFile;

fprintf('Splitting chiasmata...\n');
out = connectEM.splitChiasmataMultiSuper(param, axonFile, {splitFile});
out.info = info;

%% saving result
Util.saveStruct(outFile, out);
system(sprintf('chmod a-w "%s"', outFile));

%% utilities
function splitData = buildSplitData(taskGenFile, taskIdFile, nmlDir)
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
end