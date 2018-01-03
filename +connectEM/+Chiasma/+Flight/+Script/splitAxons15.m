% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.Flight.loadSplitData;
import connectEM.Chiasma.Flight.shuffleExits;

clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171219T214041-on-axons-15a');

% task specifics
taskGenIds = {'20180102T184038'; '20171219T221304'};
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

%% running chiasma splitting
axonFile = load(tasks(end).splitDataFile, 'axonFile');
axonFile = axonFile.axonFile;

chiasmataFile = load(tasks(end).genFile, 'info');
chiasmataFile = chiasmataFile.info.param.chiasmataFile;

fprintf('Splitting chiasmata...\n');
[splitAxons, openExits] = ...
    connectEM.splitChiasmataMultiSuper( ...
        param, chiasmataFile, axonFile, {tasks.splitDataFile});
splitAxons.info = info;

%% saving result
Util.saveStruct(splitAxonsFile, splitAxons);
system(sprintf('chmod a-w "%s"', splitAxonsFile));

%% requerying
taskGen = load(tasks(end).genFile);

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
    ismember(taskGen.exits, openExits, 'rows') | ismember( ...
    taskGen.exits(:, {'aggloId', 'chiasmaId'}), unsolvedT, 'rows');
clear unsolvedT;

requery = struct;
requery.exits = taskGen.exits(requeryMask, :);
requery.taskDefs = taskGen.taskDefs(requeryMask, :);
clear requeryMask;

[requery.exits, shuffledRows] = shuffleExits(requery.exits, 500);
requery.taskDefs = requery.taskDefs(shuffledRows, :);
clear shuffledRows;

requeryTaskDefFile = fullfile( ...
    taskGenDir, sprintf('%s_flightTasks.txt', runId));
writetable( ...
    requery.taskDefs, requeryTaskDefFile, 'WriteVariableNames', false);

%% build output
requery.info = info;
requery.taskDefFile = requeryTaskDefFile;
requery.chiasmataFile = taskGen.chiasmataFile;

requeryTaskGenFile = fullfile( ...
    taskGenDir, sprintf('%s_taskGeneration.mat', runId));
Util.saveStruct(requeryTaskGenFile, requery);
