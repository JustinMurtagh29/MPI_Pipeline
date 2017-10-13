% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

clear;
outputDir = '/home/amotta/Desktop';
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
load(fullfile(rootDir, 'allParameter.mat'), 'p');

oldAxons = load(fullfile( ...
    p.saveFolder, 'aggloState', 'axons_06_c.mat'));
axons = oldAxons.axons(oldAxons.indBigAxons);

curDir = fullfile( ...
    p.saveFolder, 'chiasmataSplitting', ...
    '20171009T193744-kmb-on-axons-6c');

curData = load(fullfile(curDir, 'input-data.mat'));

queries = table;
queries.axonId = curData.queries(:, 1);
queries.chiasmaId = curData.queries(:, 2);
queries.exitId = curData.queries(:, 3);
queries.seedPos = curData.queries(:, 4:6);
queries.centerNodeId = curData.queries(:, 7);
queries.taskId = curData.taskIds;
flights = curData.ff;

clear curDir curData;

fprintf('# handed out queries: %d\n', size(queries, 1));
fprintf('# queries answered: %d\n', numel(flights.segIds));

%%
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

% Clean-up
queries(:, {'flightId', 'uniChiasmaId'}) = [];

%% Limit to 4-fold chiasmata
[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');

uniChiasma4Fold = find(accumarray(queries.uniChiasmaId, 1) == 4);
fprintf('# 4-fold chiasmata answered: %d\n', numel(uniChiasma4Fold));

queries = queries(ismember( ...
    queries.uniChiasmaId, uniChiasma4Fold), :);
queries.uniChiasmaId = [];

%%
% Group queries by axons
queries.row = reshape(1:size(queries, 1), [], 1);
[uniAxonIds, ~, queries.uniAxonId] = unique(queries.axonId);

% Start applying chiasmata
queries.overlaps(:) = {[]};

for curAxonIdx = 1:numel(uniAxonIds)
    curAxon = axons(uniAxonIds(curAxonIdx));
    curAxon.nodesScaled = bsxfun(@times, ...
        curAxon.nodes(:, 1:3), p.raw.voxelSize);
    curQueries = queries( ...
        queries.uniAxonId == curAxonIdx, :);
    
   [~, curSummary] = ...
        connectEM.splitChiasmataMulti( ...
            p, curAxon, curQueries, 'dryRun', true);
        
    curOverlaps = cellfun(@(t) ...
        t.overlaps, curSummary.tracings, ...
        'UniformOutput', false);
    curOverlaps = cat(1, curOverlaps{:});
    
    queries.overlaps(curQueries.row) = curOverlaps;
end

queries.uniAxonId = [];

%% evaluate flight path attachment
numEndings = cellfun(@numel, queries.overlaps);

flightEval = table;
flightEval.twoEndings = (numEndings == 2);
flightEval.tooFewEndings = (numEndings < 2);
flightEval.tooManyEndings = (numEndings > 2);

temp = varfun(@sum, flightEval);
temp.Properties.VariableNames = flightEval.Properties.VariableNames;

fprintf('\n');
fprintf('Flight path evaluation:\n\n');
disp(temp);

%% Look at flight paths with too few attachments
% Make sure that start was correctly identified
tooFewQueries = queries(flightEval.tooFewEndings, :);
mean(arrayfun( ...
    @(a, b) ismember(a, b{1}), ...
    tooFewQueries.exitId, tooFewQueries.overlaps));

rng(0);
randIds = randperm(size(tooFewQueries, 1));
randIds = randIds(1:20);

curOutputDir = fullfile(outputDir, 'too-few-endings');
if ~exist(curOutputDir, 'dir'); mkdir(curOutputDir); end

for curIdx = randIds
    skel = skeleton();
    
    curQuery = tooFewQueries(curIdx, :);
    curName = curQuery.taskId;
    curNodes = reshape(curQuery.flightNodes{1}, [], 3);
    skel = skel.addTree(curName, curNodes);
    
    % show super-agglomerate
    curAxon = axons(curQuery.axonId);
    skel = skel.addTree('Axon', curAxon.nodes(:, 1:3), curAxon.edges);
    
    % show other exits of chiasma
    curQueries = queries( ...
        queries.axonId == curQuery.axonId ...
      & queries.chiasmaId == curQuery.chiasmaId, :);
    curCenterPos = curAxon.nodes(curQuery.centerNodeId, 1:3);
    
    for curExitIdx = 1:size(curQueries, 1)
        skel = skel.addTree( ...
            sprintf('Exit %d', curQueries.exitId(curExitIdx)), ...
            cat(1, curQueries.seedPos(curExitIdx, :), curCenterPos));
    end
    
    skel = Skeleton.setParams4Pipeline(skel, p);
    skel.write(fullfile(curOutputDir, ...
        sprintf('too-few-endings-%d.nml', curIdx)));
end

%%
[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');
chiasmaCount = max(queries.uniChiasmaId);
chiasmaEval = array2table( ...
    false(chiasmaCount, 3), 'VariableNames', ...
    {'invalidOverlaps', 'invalidBidir', 'valid'});

for curIdx = 1:chiasmaCount
    curQueries = queries(queries.uniChiasmaId == curIdx, :);
    curOverlaps = curQueries.overlaps;
    
    if ~all(cellfun(@numel, curOverlaps) == 2)
        % all flight paths must have start and end
        chiasmaEval.invalidOverlaps(curIdx) = true;
        continue;
    end
    
    curOverlaps = cell2mat(curOverlaps(:)');
    curOverlaps = sort(curOverlaps, 1)';
   [curPairs, ~, curPairCount] = unique(curOverlaps, 'rows');
    curPairCount = accumarray(curPairCount, 1);
    
    if ~all(curPairCount == 2)
        % flight paths must be in mutual agreement
        chiasmaEval.invalidBidir(curIdx) = true;
        continue;
    end
    
    chiasmaEval.valid(curIdx) = true;
end

temp = varfun(@sum, chiasmaEval);
temp.Properties.VariableNames = chiasmaEval.Properties.VariableNames;

fprintf('4-fold chiasmata evaluation:\n\n');
disp(temp);

%% show chiasmata which invalidate bidirectionality
chiasmaIds = find(chiasmaEval.invalidBidir);
chiasmaIds = reshape(chiasmaIds, 1, []);

rng(0);
randIds = randperm(numel(chiasmaIds));
randIds = chiasmaIds(randIds(1:10));

curOutputDir = fullfile(outputDir, 'invalid-bidir');
if ~exist(curOutputDir, 'dir'); mkdir(curOutputDir); end

for curIdx = randIds
    skel = skeleton();
    
    curQueries = queries(queries.uniChiasmaId == curIdx, :);
    curAxon = axons(curQueries.axonId(1));
    
    skel = skel.addTree('Axon', curAxon.nodes(:, 1:3), curAxon.edges);
    
    for curQueryIdx = 1:size(curQueries, 1)
        curExitId = curQueries.exitId(curQueryIdx);
        curOverlaps = curQueries.overlaps{curQueryIdx};
        
        assert(any(curOverlaps == curExitId));
        curOverlaps(curOverlaps == curExitId) = [];
        
        curQueryName = sprintf( ...
            'Flight %d → %d', curExitId, curOverlaps(1));
        skel = skel.addTree( ...
            curQueryName, curQueries.flightNodes{curQueryIdx});
    end
    
    curCenterPos = curAxon.nodes(curQueries.centerNodeId(1), 1:3);
    
    for curExitIdx = 1:size(curQueries, 1)
        skel = skel.addTree( ...
            sprintf('Exit %d', curQueries.exitId(curExitIdx)), ...
            cat(1, curQueries.seedPos(curExitIdx, :), curCenterPos));
    end
    
    skel = Skeleton.setParams4Pipeline(skel, p);
    skel.write(fullfile(curOutputDir, ...
        sprintf('invalid-bidir-%d.nml', curIdx)));
end

%% calculate partitioning
[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');
chiasmaCount = max(queries.uniChiasmaId);
chiasmaPartition = cell(chiasmaCount, 1);

for curIdx = 1:chiasmaCount
    curQueries = queries(queries.uniChiasmaId == curIdx, :);
    curNrExits = size(curQueries, 1);
    
    curOverlaps = curQueries.overlaps;
    if ~all(cellfun(@numel, curOverlaps) <= 2)
        % Do not consider chiasmata which contain flight paths that attach
        % to too many exits. Having less than two exit overlaps is okay.
        % This may happen in 3 + 1, 2 + 1 + 1 or 1 + 1 + 1 + 1 partitions.
        chiasmaPartition{curIdx} = zeros(curNrExits, 1);
        continue;
    end
    
    % Pad overlaps with zeros
    % → Non-attachment is exit #0
    curOverlaps = cellfun( ...
        @(o) cat(1, o, zeros(2 - numel(o), 1)), ...
        curOverlaps, 'UniformOutput', false);
    curPairsAll = sort(cat(2, curOverlaps{:}), 1)';
    
    curPairs = curPairsAll;
    curPairs(~all(curPairs, 2), :) = [];
    curPairs = unique(curPairs, 'rows');
    curPairs = reshape(curPairs, [], 2);
    
    % Find parition
    curAdj = sparse( ...
        curPairs(:, 2), curPairs(:, 1), ...
        true, curNrExits, curNrExits);
    
   [~, curLUT] = graphconncomp(curAdj, 'Directed', false);
   
    % Build partition
    curPartition = accumarray(curLUT(:), 1);
    curPartition = sort(curPartition, 'descend');
    
    chiasmaPartition{curIdx} = curPartition;
end

%% evaluate partitions
chiasmaPartitionStr = cellfun(@(p) ...
    strjoin(arrayfun(@num2str, p(:)', 'Uni', false), '-'), ...
    chiasmaPartition, 'UniformOutput', false);
[uniPartitions, ~, uniPartitionCount] = unique(chiasmaPartitionStr);
uniPartitionCount = accumarray(uniPartitionCount, 1);


partitionEval = table;
partitionEval.partition = uniPartitions;
partitionEval.count = uniPartitionCount;

fprintf('Evaluation of chiasmata partitioning:\n\n');
disp(partitionEval);

%% find 3 + 1 chiasmata
curPartitionStr = '3-1';
chiasmaIds = find(strcmpi(chiasmaPartitionStr, curPartitionStr));
chiasmaIds = reshape(chiasmaIds, 1, []);

rng(0);
randIds = randperm(numel(chiasmaIds));
randIds = chiasmaIds(randIds(1:min(10, numel(chiasmaIds))));

curOutputDir = fullfile( ...
    outputDir, sprintf('%s-partition', curPartitionStr));
if ~exist(curOutputDir, 'dir'); mkdir(curOutputDir); end

for curIdx = randIds
    skel = skeleton();
    
    curQueries = queries(queries.uniChiasmaId == curIdx, :);
    curAxon = axons(curQueries.axonId(1));
    
    skel = skel.addTree('Axon', curAxon.nodes(:, 1:3), curAxon.edges);
    
    for curQueryIdx = 1:size(curQueries, 1)
        curExitId = curQueries.exitId(curQueryIdx);
        curOverlaps = curQueries.overlaps{curQueryIdx};
        
        assert(any(curOverlaps == curExitId));
        curOverlaps(curOverlaps == curExitId) = [];
        curOverlaps(end + 1) = 0; %#ok
        
        curQueryName = sprintf( ...
            'Flight %d → %d', curExitId, curOverlaps(1));
        skel = skel.addTree( ...
            curQueryName, curQueries.flightNodes{curQueryIdx});
    end
    
    curCenterPos = curAxon.nodes(curQueries.centerNodeId(1), 1:3);
    
    for curExitIdx = 1:size(curQueries, 1)
        skel = skel.addTree( ...
            sprintf('Exit %d', curQueries.exitId(curExitIdx)), ...
            cat(1, curQueries.seedPos(curExitIdx, :), curCenterPos));
    end
    
    skel = Skeleton.setParams4Pipeline(skel, p);
    skel.write(fullfile(curOutputDir, sprintf( ...
        '%s-partition-%d.nml', curPartitionStr, curIdx)));
end