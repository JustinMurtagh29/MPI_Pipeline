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
axonFile = tasks.info.param.axonFile;
chiParam = tasks.info.param.chiParam;
tasks = struct2table(tasks.taskDefs);

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

axons  = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons(:));
axons = axons.axons(axonIds);

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
[queries.exits, queries.comments] = cellfun( ...
    @connectEM.Chiasma.Ortho.parseQuery, ...
    queries.nmlFile, 'UniformOutput', false);

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
commentMask = ~cellfun(@isempty, queries.comments);
missingExitMask = ~cellfun( ...
    @(a, b) size(a, 1) == size(b, 1), ...
    queries.exits, queries.exitNodeIds);

fprintf('\n');
fprintf('# queries answered: %d\n', size(queries, 1));
fprintf('# queries with comment: %d\n', sum(commentMask));
fprintf('# queries with missing exits: %d\n', sum(missingExitMask));

fprintf('⇒ Removing these queries...\n');
queries(commentMask | missingExitMask, :) = [];
clear commentMask missingExitMask;

%% apply results