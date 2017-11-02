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

info = Util.runInfo();

%% load tasks and task definitions
tasks = load(taskFile);

% load axons
axonFile = tasks.info.param.axonFile;
axons  = load(axonFile, 'axons', 'indBigAxons');

% restrict to large axons
axonIds = find(axons.indBigAxons(:));
axons = axons.axons(axonIds);

% get parameters for chiasma splitting
chiParam = tasks.info.param.chiParam;
chiParam.sphereRadiusOuter = inf;

% get task definitions
tasks = struct2table(tasks.taskDefs);

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

tasksIds = connectEM.Chiasma.Ortho.loadTaskIds(taskIdFile);

%% replace `nmlFile` by `taskId`
[~, tasks.id] = ismember(tasks.nmlFile, tasksIds.nmlFile);
tasks.id = tasksIds.id(tasks.id);

tasks.nmlFile = [];
tasks = circshift(tasks, 1, 2);

%% find NML files
nmlFiles = NML.findFiles(nmlDir);
nmlFiles = reshape(nmlFiles, [], 1);

queries = table;
queries.nmlFile = fullfile(nmlDir, nmlFiles);

% parse NML files
[queries.exits, queries.valid] = cellfun( ...
    @connectEM.Chiasma.Ortho.parseQuery, ...
    queries.nmlFile, 'UniformOutput', false);

queriesValid = cat(1, queries.valid{:});
queriesValid = struct2table(queriesValid);

queries.valid = [];
queries = cat(2, queries, queriesValid);
clear queriesValid;

%% add task data to queries
% extract task IDs
queryData = cellfun( ...
    @(s) strsplit(s, '__'), ...
    nmlFiles, 'UniformOutput', false);
queryData = cat(1, queryData{:});
queries.taskId = queryData(:, 2);
clear queryData;

[~, taskRow] = ismember(queries.taskId, tasks.id);
queries = cat(2, queries, tasks(taskRow, :));
queries(:, {'id', 'taskId'}) = [];
clear taskRow;

%% statistics
treeMask = queries.nrInvalidTrees > 0;
commentMask = queries.nrInvalidComments > 0;
exitMismatchMask = cellfun( ...
    @(a, b) size(a, 1) ~= size(b, 1), ...
    queries.exits, queries.exitNodeIds);

fprintf('\n');
fprintf('# queries answered: %d\n', size(queries, 1));
fprintf('# queries with invalid tree: %d\n', sum(treeMask));
fprintf('# queries with invalid comment: %d\n', sum(commentMask));
fprintf('# queries with invalid exits: %d\n', sum(exitMismatchMask));

fprintf('⇒ Removing these queries...\n');
queries(treeMask | commentMask | exitMismatchMask, :) = [];
clear commentMask missingExitMask;

fprintf('\n');
fprintf('# queries left: %d\n', size(queries, 1));

% sort by axon
queries = sortrows(queries, 'axonId');

%% apply results
uniAxonIds = unique(queries.axonId);

for curIdx = 1:numel(uniAxonIds)
    curAxonId = uniAxonIds(curIdx);
    curMask = queries.axonId == curAxonId;
    
    curAxons = connectEM.Chiasma.Ortho.splitWithQueries( ...
        param, chiParam, axons(curAxonId), queries(curMask, :));
end