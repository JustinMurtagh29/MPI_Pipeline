% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

inputVersion = '20171101T173904';
inputDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171031T000000-ortho-mode-on-axons-8');

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

% load pipeline parameters
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% load task IDs and queries
taskIds = connectEM.Chiasma.Util.loadTaskIds(taskIdFile);
queries = connectEM.Chiasma.Ortho.loadQueries(nmlDir);

%% split chiasmata
out = connectEM.Chiasma.Ortho.split( ...
    param, chiParam, oldAxons, tasks, taskIds, queries);

%% save output
out.info = info;

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

outFile = sprintf('%s_results.mat', datestr(now, 30));
outFile = fullfile(outputDir, outFile);
Util.saveStruct(outFile, out);
