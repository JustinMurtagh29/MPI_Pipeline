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

chiasmaDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/round4';

% task specifics
taskGenIds = {'20171115T135555'};
taskGenDir = fullfile(chiasmaDir, 'taskGeneration');

tasks = cellfun(@(id) struct( ...
    'genFile', fullfile( ...
        taskGenDir, sprintf('%s_taskGeneration.mat', id)), ...
    'idFile', fullfile( ...
        chiasmaDir, 'taskAnswers', sprintf('%s_flightTaskIDs.txt', id)), ...
    'nmlDir', fullfile( ...
        chiasmaDir, 'taskAnswers', sprintf('%s_flightTasks', id)), ...
    'splitDataFile', fullfile( ...
        chiasmaDir, 'chiasmataSplitting', sprintf('%s_splitData.mat', id)), ...
	'splitDataFileBuild', false), ...
    taskGenIds);

% output file
splitAxonsFile = fullfile( ...
    chiasmaDir, sprintf('%s_splitAxons.mat', taskGenIds{1}));

chiParam = load('/tmpscratch/scchr/L4/chiasmata/dendrites/round4/round4_inclInfos_chiasmata.mat','info');
chiParam = chiParam.info.param.chiasmaParam;
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

fprintf('Splitting chiasmata...\n');
[splitAxons, openExits] = ...
    connectEM.splitChiasmataMultiSuper( ...
        param, chiParam, axonFile, {tasks.splitDataFile}, 1);
splitAxons.info = info;

%% saving result
Util.saveStruct(splitAxonsFile, splitAxons);
system(sprintf('chmod a-w "%s"', splitAxonsFile));

%% requerying
taskGen = load(tasks(end).genFile);
exits = taskGen.exits;
taskDefs = taskGen.taskDefs;
chiasmataFile = taskGen.chiasmataFile;
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
requery.taskDefs = taskDefs;

% NOTE(amotta): Martin and I independently discovered a huge bug in the
% generation of requeries. This bug was present throughout the entire L4
% project, and rendered requeries effectively useless.
%   To properly document the L4 project, the bug is left in place. But we
% prevent anybody from using this code without uncommenting the fix first.
error('Uncomment the following bug-fix before generating requeries');
% requery.taskDefs = taskDefs(requeryMask, :);
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
requery.chiasmataFile = chiasmataFile;

requeryTaskGenFile = fullfile( ...
    taskGenDir, sprintf('%s_taskGeneration.mat', runId));
Util.saveStruct(requeryTaskGenFile, requery);
