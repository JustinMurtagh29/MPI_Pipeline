% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.*;
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
taskGenIds = {'20171201T144024'; '20171128T155724'; '20171124T161031'};
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

%% determine overlaps
fprintf('Splitting chiasmata...\n');
[dryRun, openExits] = ...
    connectEM.splitChiasmataMultiSuper( ...
        param, chiasmataFile, ...
        axonFile, {tasks.splitDataFile}, ...
        'dryRun', true, 'partialAnswers', true);

axons = dryRun.oldAxons.axons;
bigAxonMask = dryRun.oldAxons.indBigAxons;
bigAxonIds = find(bigAxonMask);

% fix empty edges
newEdges = arrayfun( ...
    @(a) reshape(a.edges, [], 2), ...
    axons, 'UniformOutput', false);
[axons.edges] = deal(newEdges{:});
clear newEdges;

bigAxons = axons(bigAxonIds);

%% do auto-completion of overlaps
chiasmata = load(chiasmataFile, 'info', 'chiasmata');
chiasmaParam = chiasmata.info.param.chiasmaParam;
chiasmata = chiasmata.chiasmata;

overlaps = Triplet.buildOverlaps(chiasmata, dryRun.summary);

%% split axons
fprintf('Splitting agglomerates... ');
splitAxons = arrayfun(@(a, c, o) ...
    Triplet.splitAgglo(param, chiasmaParam, a, c{1}, o{1}), ...
    bigAxons, chiasmata, overlaps, 'UniformOutput', false);
fprintf('done!\n');

out = struct;
out.axons = cat(1, splitAxons{:});
out.parentIds = repelem(bigAxonIds, cellfun(@numel, splitAxons));

% add small agglomerates
otherAxonIds = setdiff((1:numel(axons))', bigAxonIds);
out.axons = cat(1, out.axons, axons(otherAxonIds));
out.parentIds = cat(1, out.parentIds, otherAxonIds);

% build `indBigAxons` mask
out.indBigAxons = bigAxonMask(out.parentIds);

% add info
out.info = info;

%% saving result
Util.saveStruct(splitAxonsFile, out);
system(sprintf('chmod a-w "%s"', splitAxonsFile));

%% debugging
%{
flights = arrayfun( ...
    @(t) subsref( ...
        load(t.splitDataFile), ...
        struct('type', '.', 'subs', 'ff')), ...
    tasks, 'UniformOutput', false);
flights = Util.concatStructs(1, flights{:});

rng(0);
debugAxonIds = cat(1, dryRun.summary.axonId);
debugAxonIds = debugAxonIds(randperm(numel(debugAxonIds), 10));

debugSkels = Triplet.debug( ...
    bigAxons, splitAxons, dryRun.summary, flights, debugAxonIds);

for curIdx = 1:numel(debugAxonIds)
    curSkel = debugSkels{curIdx};
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curSkelName = sprintf('axon-%d.nml', debugAxonIds(curIdx));
    curSkel.write(fullfile('/home/amotta/Desktop/debug', curSkelName));
end
%}

%% generate next round of queries
% restrict to unsolved triplets
%{
chiasmaT = Detect.buildTable(chiasmata, bigAxons);
chiasmaT(chiasmaT.nrExits ~= 3, :) = [];
chiasmaT(chiasmaT.isSolved, :) = [];

% select exits for next query round
exits = Flight.selectExits(chiasmaT, overlaps, inf);

taskDefFile = fullfile(taskGenDir, sprintf('%s_flightTasks.txt', runId));
taskDefs = Flight.generateTasks(param, chiasmata, exits, taskDefFile);

%% build output
out = struct;
out.info = info;
out.exits = exits;
out.taskDefs = taskDefs;
out.taskDefFile = taskDefFile;
out.chiasmataFile = chiasmataFile;

outFile = sprintf('%s_taskGeneration.mat', runId);
Util.saveStruct(fullfile(taskGenDir, outFile), out);
%}