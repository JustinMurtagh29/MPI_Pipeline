% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.Flight.loadSplitData;
import connectEM.Chiasma.Flight.shuffleExits;
import connectEM.Chiasma.Util.loadFlightPaths;
import connectEM.Chiasma.Util.loadTaskIds;

clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a');

% task specifics
taskGenId = '20171108T105021';
taskGenDir = fullfile(chiasmaDir, 'taskGeneration');

taskGenFile = fullfile( ...
    taskGenDir, ...
    sprintf('%s_taskGeneration.mat', taskGenId));
taskIdFile = fullfile( ...
    chiasmaDir, 'taskAnswers', ...
    sprintf('%s_flightTaskIDs.txt', taskGenId));
nmlDir = fullfile( ...
    chiasmaDir, 'taskAnswers', ...
    sprintf('%s_flightTasks', taskGenId));

splitDataFile = fullfile( ...
    chiasmaDir, sprintf('%s_splitData.mat', taskGenId));
splitDataFileBuild = false;

% output file
splitAxonsFile = fullfile( ...
    chiasmaDir, sprintf('%s_splitAxons.mat', taskGenId));

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% build or load split file
if splitDataFileBuild || ~exist(splitDataFile, 'file')
    splitData = loadSplitData(taskGenFile, taskIdFile, nmlDir);
    splitData.info = info;
    
    Util.saveStruct(splitDataFile, splitData);
    clear splitData;
end

%% running chiasma splitting
axonFile = load(splitDataFile, 'axonFile');
axonFile = axonFile.axonFile;

fprintf('Splitting chiasmata...\n');
[splitAxons, openExits] = ...
    connectEM.splitChiasmataMultiSuper(param, axonFile, {splitDataFile});
splitAxons.info = info;

%% saving result
Util.saveStruct(splitAxonsFile, splitAxons);
system(sprintf('chmod a-w "%s"', splitAxonsFile));

%% requerying
taskGen = load(taskGenFile);
exits = taskGen.exits;
taskDefs = taskGen.taskDefs;
clear taskGen;

% find unsolved chiasmata
unsolvedT = table;
unsolvedT.aggloId = repelem( ...
    cat(1, splitAxons.summary.axonId), ...
    cat(1, splitAxons.summary.nrChiasmata));
unsolvedT.chiasmaId = cat(1, splitAxons.summary.chiasmaId);

unsolvedT.isSolved = cat(1, splitAxons.summary.solved);
unsolvedT(unsolvedT.isSolved, :) = [];
unsolvedT.isSolved = [];

% mark tasks to requery
requeryMask =  ...
    ismember(exits, openExits, 'rows') | ismember( ...
    exits(:, {'aggloId', 'chiasmaId'}), unsolvedT, 'rows');
clear unsolvedT;

requery = struct;
requery.exits = exits(requeryMask, :);
clear requeryMask;

[requery.exits, shuffledRows] = shuffleExits(requery.exits, 500);
requery.taskDefs = taskDefs(shuffledRows, :);
clear shuffledRows;

requeryTaskDefFile = fullfile( ...
    taskGenDir, sprintf('%s_flightTasks.txt', runId));
writetable( ...
    requery.taskDefs, requeryTaskDefFile, 'WriteVariableNames', false);

%% build output
requery.info = info;
requery.taskDefFile = requeryTaskDefFile;

requeryTaskGenFile = fullfile( ...
    taskGenDir, sprintf('%s_taskGeneration.mat', runId));
Util.saveStruct(requeryTaskGenFile, requery);