% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.Flight.selectExits;
import connectEM.Chiasma.Flight.generateTasks;

%%
clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a');

chiasmataFile = fullfile(chiasmaDir, '20171104T184018_chiasmata.mat');
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

%% build empty overlaps
overlaps = cellfun(@(c) cellfun(@(q) ...
    zeros(numel(q), 2), c.queryIdx, 'UniformOutput', false), ...
	chiasmata, 'UniformOutput', false);

%% find exits to query
exits = selectExits(axons, chiasmata, overlaps);

taskDefFile = fullfile(outputDir, sprintf('%s_flightTasks.txt', runId));
taskDefs = generateTasks(param, chiasmata, exits, taskDefFile);

%% build output
out = struct;
out.info = info;
out.exits = exits;
out.taskDefs = taskDefs;
out.taskDefFile = taskDefFile;

outFile = sprintf('%s_taskGeneration.mat', runId);
Util.saveStruct(fullfile(outputDir, outFile), out);