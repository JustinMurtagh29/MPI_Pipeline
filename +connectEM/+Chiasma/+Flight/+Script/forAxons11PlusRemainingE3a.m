% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.Detect.buildTable;
import connectEM.Chiasma.Flight.selectExits;
import connectEM.Chiasma.Flight.shuffleExits;
import connectEM.Chiasma.Flight.generateTasks;

%%
clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171113T171102-on-axons-11a-plus-remaining-E3a');

chiasmataFile = fullfile(chiasmaDir, '20171117T174027_chiasmata.mat');
outputDir = fullfile(chiasmaDir, 'taskGeneration');
clear workingDir;

info = Util.runInfo();

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

chiasmata = load(chiasmataFile);
axonFile = chiasmata.info.param.axonFile;
chiasmata = chiasmata.chiasmata;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% find exits to query
chiasmaT = buildTable(chiasmata, axons);
chiasmaT(chiasmaT.isSolved, :) = [];

exits = selectExits(chiasmaT, [], inf);
exits = shuffleExits(exits, 500);

taskDefFile = fullfile(outputDir, sprintf('%s_flightTasks.txt', runId));
taskDefs = generateTasks(param, chiasmata, exits, taskDefFile);

%% build output
out = struct;
out.info = info;
out.exits = exits;
out.taskDefs = taskDefs;
out.taskDefFile = taskDefFile;
out.chiasmataFile = chiasmataFile;

outFile = sprintf('%s_taskGeneration.mat', runId);
Util.saveStruct(fullfile(outputDir, outFile), out);