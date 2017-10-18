% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

clear;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
load(fullfile(rootDir, 'allParameter.mat'), 'p');

oldAxons = load(fullfile( ...
    p.saveFolder, 'aggloState', 'axons_06_c.mat'));

bigAxonIds = find(oldAxons.indBigAxons(:));
axons = oldAxons.axons(bigAxonIds);

curDir = fullfile( ...
    p.saveFolder, 'chiasmataSplitting', ...
    '20171009T193744-kmb-on-axons-6c');
outputDir = fullfile(curDir, 'outputs');
curData = load(fullfile( ...
    curDir, '20171018T104038_input-data.mat'));

queries = table;
queries.axonId = curData.queries(:, 1);
queries.chiasmaId = curData.queries(:, 2);
queries.exitId = curData.queries(:, 3);
queries.exitNodeId = curData.queries(:, 8);
queries.seedPos = curData.queries(:, 4:6);
queries.centerNodeId = curData.queries(:, 7);
queries.taskId = curData.taskIds;
flights = curData.ff;

clear curDir curData;

fprintf('# handed out queries: %d\n', size(queries, 1));
fprintf('# queries answered: %d\n', numel(flights.segIds));

%% cleaning up input data
% Sort queries to make processing easier
queries = sortrows(queries, {'axonId', 'chiasmaId', 'exitId'});
flights = structfun(@(x) x(:), flights, 'UniformOutput', false);

% Assign flight paths to queries
[~, queries.flightId] = ismember( ...
    queries.taskId, flights.filenamesShort);

% Find chiasmata for which we have all queries answered
[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');
uniChiasmaDoneIds = find(accumarray( ...
    queries.uniChiasmaId, queries.flightId, [], @all));
fprintf('# chiasmata answered: %d\n', numel(uniChiasmaDoneIds));

% Limit ourselves to done chiasmata
queries = queries(ismember( ...
    queries.uniChiasmaId, uniChiasmaDoneIds), :);

queries.flightNodes = flights.nodes(queries.flightId);
queries.flightSegIds = flights.segIds(queries.flightId);
queries.flightComment = flights.comments(queries.flightId);

fprintf('\n');
fprintf('# queries without nodes: %d\n', ...
    sum(cellfun(@isempty, queries.flightNodes)));
fprintf('# queries with comment: %d\n', ...
    sum(~cellfun(@isempty, queries.flightComment)));
fprintf('⇒ Removing the above queries ...\n');

% removing invalid flights and their chiasmata
invalidFlightMask = ...
    cellfun(@isempty, queries.flightNodes) ...
 | ~cellfun(@isempty, queries.flightComment);

invalidChiasmaIds = unique(queries.uniChiasmaId(invalidFlightMask));
queries(ismember(queries.uniChiasmaId, invalidChiasmaIds), :) = [];

[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');
uniChiasmaDoneIds = find(accumarray( ...
    queries.uniChiasmaId, queries.flightId, [], @all));

fprintf('\n');
fprintf('# answered chiasmata left: %d\n', numel(uniChiasmaDoneIds));

% Clean-up
queries(:, {'flightId', 'uniChiasmaId'}) = [];

%% process axons
% Group queries by axons
queries.execute = false(size(queries, 1), 1);
queries.row = reshape(1:size(queries, 1), [], 1);
[uniAxonIds, ~, queries.uniAxonId] = unique(queries.axonId);

% Start applying chiasmata
queries.overlaps(:) = {[]};

axonsSplit = cell(numel(uniAxonIds), 1);
summaries = cell(numel(uniAxonIds), 1);

tic
for curAxonIdx = 1:numel(uniAxonIds)
    curAxon = axons(uniAxonIds(curAxonIdx));
    curAxon.nodesScaled = bsxfun(@times, ...
        curAxon.nodes(:, 1:3), p.raw.voxelSize);
    curQueries = queries( ...
        queries.uniAxonId == curAxonIdx, :);
    
   [axonsSplit{curAxonIdx}, curSummary] = ...
        connectEM.splitChiasmataMulti(p, curAxon, curQueries);
        
    curOverlaps = cellfun(@(t) ...
        t.overlaps, curSummary.tracings, ...
        'UniformOutput', false);
    curOverlaps = cat(1, curOverlaps{:});
    
    summaries{curAxonIdx} = curSummary;
    queries.overlaps(curQueries.row) = curOverlaps;
end
toc

summaries = cat(1, summaries{:});
queries.uniAxonId = [];

%% evaluate flight path attachment
numEndings = cellfun(@(o) sum(o > 0), queries.overlaps);

flightEval = table;
flightEval.twoEndings = (numEndings == 2);
flightEval.tooFewEndings = (numEndings < 2);

temp = varfun(@sum, flightEval);
temp.Properties.VariableNames = flightEval.Properties.VariableNames;

fprintf('\n');
fprintf('Flight path evaluation:\n\n');
disp(temp);

%% evaluate partitions
chiasma = arrayfun(@(s) ...
    cellfun(@(t) { ...
        t.partition(:), t.partitionIsValid}, ...
        s.tracings, 'UniformOutput', false), ...
    summaries, 'UniformOutput', false);
chiasma = cat(1, chiasma{:});
chiasma = cat(1, chiasma{:});

chiasma = cell2table( ...
    chiasma, 'VariableNames', {'partition', 'valid'});

makePartitionStr = @(p) strjoin( ...
    arrayfun(@num2str, p(:)', 'Uni', false), '-');
chiasmaPartitionStr = cellfun( ...
    makePartitionStr, chiasma.partition, 'UniformOutput', false);
[uniPartitions, ~, uniPartitionIds] = unique(chiasmaPartitionStr);

partitionEval = table;
partitionEval.partition = uniPartitions;
partitionEval.count = accumarray(uniPartitionIds, 1);
partitionEval.valid = accumarray(uniPartitionIds, chiasma.valid);

fprintf('Evaluation of chiasmata partitioning:\n\n');
disp(partitionEval);

%%
% NOTE(amotta): Add `endings` field to old agglomerates
[oldAxons.axons.endings] = deal([]);
oldAxons.indBigAxons = oldAxons.indBigAxons(:);

out = struct;
out.p = p;
out.oldAxons = oldAxons;
out.gitInfo = Util.gitInfo();

out.summary = summaries;
out.summaryIds = bigAxonIds(uniAxonIds);

out.axons = cat(1, axonsSplit{:});
out.parentIds = repelem(out.summaryIds, cellfun(@numel, axonsSplit));
otherAxons = setdiff((1:numel(oldAxons.axons))', out.summaryIds);

% add small agglomerates
out.axons = cat(1, out.axons, oldAxons.axons(otherAxons));
out.parentIds = cat(1, out.parentIds, otherAxons);

% build `indBigAxons` mask
out.indBigAxons = oldAxons.indBigAxons(out.parentIds);

outFile = sprintf('%s_results.mat', datestr(now, 30));
outFile = fullfile(outputDir, outFile);
Util.saveStruct(outFile, out);
