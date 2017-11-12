% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
import connectEM.Chiasma.Flight.loadSplitData;
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
