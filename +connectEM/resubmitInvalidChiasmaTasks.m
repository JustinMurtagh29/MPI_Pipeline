% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

clear;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

dataDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171009T193744-kmb-on-axons-6c');
outputDir = fullfile( ...
    dataDir, 'requeries', 'raw-data');
data = load(fullfile( ...
    dataDir, '20171018T205736_input-data.mat'));

queries = table;
queries.axonId = data.queries(:, 1);
queries.chiasmaId = data.queries(:, 2);
queries.exitId = data.queries(:, 3);
queries.exitNodeId = data.queries(:, 8);
queries.seedPos = data.queries(:, 4:6);
queries.centerNodeId = data.queries(:, 7);
queries.taskId = data.taskIds;
queries.row = reshape( ...
    1:size(queries, 1), [], 1);
flights = data.ff;

%% find invalid tasks
% Sort queries to make processing easier
queries = sortrows(queries, {'axonId', 'chiasmaId', 'exitId'});
flights = structfun(@(x) x(:), flights, 'UniformOutput', false);

% Assign flight paths to queries
[~, queries.flightId] = ismember( ...
    queries.taskId, flights.filenamesShort);

queries.flightNodes = flights.nodes(queries.flightId);
queries.flightSegIds = flights.segIds(queries.flightId);
queries.flightComment = flights.comments(queries.flightId);

% find task without nodes or with comment
maskEmpty = cellfun(@isempty, queries.flightNodes);
maskComment = ~cellfun(@isempty, queries.flightComment);
mask = maskEmpty | maskComment;
queries = queries(mask, :);

fprintf('# queries without nodes: %d\n', sum(maskEmpty));
fprintf('# queries with comment: %d\n', sum(maskComment));

%% resubmit tasks
curDateStr = datestr(now, 30);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

out = struct;
out.gitInfo = Util.gitInfo();
out.queries = data.queries(queries.row, :);
out.taskDef = data.taskDef(queries.row, :);

saveFile = sprintf('%s_data.mat', curDateStr);
Util.saveStruct(fullfile(outputDir, saveFile), out);

taskParam = struct;
taskParam.dataSet = '2012-09-28_ex145_07x2_ROI2017';
taskParam.taskTypeId = '56d6a7c6140000d81030701e';
taskParam.expDomain = 'queriesMHchiasma';
taskParam.expMinVal = 1;
taskParam.instances = 1;
taskParam.team = 'Connectomics department';
taskParam.project = 'queriesMHchiasma';

taskDefFile = sprintf('%s_flightTasks.txt', curDateStr);
connectEM.exportTaskDefinitions( ...
    taskParam, out.taskDef, fullfile(outputDir, taskDefFile));