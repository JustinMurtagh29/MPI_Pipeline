% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Edits by
%   Christian Schramm <christian.schramm@brain.mpg.de>

import connectEM.Chiasma.Detect.buildTable;
import connectEM.Chiasma.Flight.selectExits;
import connectEM.Chiasma.Flight.shuffleExits;
import connectEM.Chiasma.Flight.generateTasks;

%%
clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/round4';
chiasmaId = '20171104T105457';
aggloState = 'dendrites_flight_02';
dendriteFile = fullfile(rootDir, 'aggloState', aggloState);
agglos = load(dendriteFile);
aggloIds = find(agglos.indBigDends);
agglos = agglos.dendrites(aggloIds);
aggloCount = numel(agglos);

chiasmataFile = fullfile(chiasmaDir, 'round4_inclInfos_chiasmata.mat');
outputDir = fullfile(chiasmaDir, 'taskGeneration');
clear workingDir;

info = Util.runInfo();

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

if exist(chiasmataFile)
    chiasmata = load(chiasmataFile);
    chiasmata = chiasmata.chiasmata;
else
    chiasmata = ...
    connectEM.Chiasma.Detect.collectResults( ...
        chiasmaDir, chiasmaId, aggloCount);
    save(chiasmataFile,'chiasmata')
end

%% find exits to query
chiasmaT = buildTable(chiasmata, agglos);
% chiasmaT(chiasmaT.isSolved, :) = [];

exits = selectExits(chiasmaT, [], inf);
exits = shuffleExits(exits, 500);

if ~exist(outputDir)
    mkdir(outputDir)
end

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