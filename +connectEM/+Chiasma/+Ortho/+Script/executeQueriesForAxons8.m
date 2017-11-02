% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
inputVersion = '20171101T173904';
inputDir = '/home/amotta/Desktop/ortho-chiasma-queries';
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

taskFile = fullfile(inputDir, sprintf('%s_tasks.mat', inputVersion));
taskIdFile = fullfile(inputDir, sprintf('%s_taskIds.txt', inputVersion));
nmlDir = fullfile(inputDir, sprintf('%s_answeredQueries', inputVersion));
outputDir = fullfile(inputDir, 'outputs');

info = Util.runInfo();

%% load tasks and task definitions
tasks = load(taskFile);

% load "old" axons
axonFile = tasks.info.param.axonFile;
oldAxons  = load(axonFile, 'axons', 'indBigAxons');

% get parameters for chiasma splitting
chiParam = tasks.info.param.chiParam;
chiParam.sphereRadiusOuter = inf;

% get task definitions
tasks = struct2table(tasks.taskDefs);

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

taskIds = connectEM.Chiasma.Ortho.loadTaskIds(taskIdFile);
queries = connectEM.Chiasma.Ortho.loadQueries(nmlDir);

%% split chiasmata
out = connectEM.Chiasma.Ortho.split( ...
    param, chiParam, oldAxons, tasks, taskIds, queries);

%% save output
outFile = sprintf('%s_results.mat', datestr(now, 30));
outFile = fullfile(outputDir, outFile);

mkdir(outputDir);
Util.saveStruct(outFile, out);
