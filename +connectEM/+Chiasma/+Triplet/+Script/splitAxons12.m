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
    rootDir, 'tripletDetection', '20171124T100638-on-axons-12a');

% task specifics
taskGenIds = {'20171124T161031'};
taskGenDir = fullfile(chiasmaDir, 'taskGeneration');

tasks = cellfun(@(id) struct( ...
    'genFile', fullfile( ...
        taskGenDir, sprintf('%s_taskGeneration.mat', id)), ...
    'idFile', fullfile( ...
        chiasmaDir, 'taskAnswers', sprintf('%s_flightTaskIDs.txt', id)), ...
    'nmlDir', fullfile( ...
        chiasmaDir, 'taskAnswers', sprintf('%s_flightTasks', id)), ...
    'splitDataFile', fullfile( ...
        chiasmaDir, sprintf('%s_splitData.mat', id)), ...
	'splitDataFileBuild', false), ...
    taskGenIds);

% output file
splitAxonsFile = fullfile( ...
    chiasmaDir, sprintf('%s_splitAxons.mat', taskGenIds{1}));

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% build or load split file
for task = reshape(tasks, 1, [])
    if task.splitDataFileBuild || ~exist(task.splitDataFile, 'file')
        splitData = loadSplitData( ...
            param, task.genFile, task.idFile, task.nmlDir);
        splitData.info = info;

        Util.saveStruct(task.splitDataFile, splitData);
        clear splitData;
    end
end

%% loading data for split
axonFile = load(tasks(end).splitDataFile, 'axonFile');
axonFile = axonFile.axonFile;

chiasmataFile = load(tasks(end).genFile, 'info');
chiasmataFile = chiasmataFile.info.param.chiasmataFile;

%%
fprintf('Splitting chiasmata...\n');
[splitAxons, openExits] = ...
    connectEM.splitChiasmataMultiSuper( ...
        param, chiasmataFile, axonFile, {tasks.splitDataFile}, ...
        'dryRun', true, 'partialAnswers', true);
