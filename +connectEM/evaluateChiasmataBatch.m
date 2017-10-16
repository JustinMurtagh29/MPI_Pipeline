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
queries.exitNodeId = curData.queries(:, 8);
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

%% Limit to 4-fold chiasmata
chiasmaFold = 4;

[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');
uniChiasmaIds = find(accumarray(queries.uniChiasmaId, 1) == chiasmaFold);

fprintf(...
    '# %d-fold chiasmata answered: %d\n', ...
    chiasmaFold, numel(uniChiasmaIds));

queries = queries(ismember( ...
    queries.uniChiasmaId, uniChiasmaIds), :);
queries.uniChiasmaId = [];

%%
%{
% Write out random flight paths. I will then manually annotate the ground
% truth overlaps. This can then be used to optimize the overlap detection
% routine.

rng(0);
randIds = randperm(size(queries, 1));
randQueries = queries(randIds(1:50), :);

curOutputDir = fullfile(outputDir, 'rand-queries');
if ~exist(curOutputDir, 'dir'); mkdir(curOutputDir); end

for curIdx = 1:size(randQueries, 1)
    skel = skeleton();
    
    curQuery = randQueries(curIdx, :);
    curNodes = reshape(curQuery.flightNodes{1}, [], 3);
    skel = skel.addTree('Flight', curNodes);
    
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
        sprintf('rand-query-%d.nml', curIdx)));
end

% extend table and manually manipulate it
randQueries.gtExit = nan(size(randQueries, 1), 1);

%%
gtFile = '/mnt/mpibr/data/Personal/mottaa/L4/2017-10-13-Chiasmata-Attachment-GT/chiasma-gt-exits.mat';
gtQueries = load(gtFile, 'randQueries');
gtQueries = gtQueries.randQueries;

% ignore entries without groundtruth annotation
gtQueries.id = reshape(1:size(gtQueries, 1), [], 1);
gtQueries(isnan(gtQueries.gtExit), :) = [];
gtCount = size(gtQueries, 1);

gtEval = table;
gtEval.id = nan(gtCount, 1);
gtEval.expected = cell(gtCount, 1);
gtEval.found = cell(gtCount, 1);
gtEval.correct = false(gtCount, 1);

for curIdx = 1:size(gtQueries, 1)
    curQuery = gtQueries(curIdx, :);
    
    curAxon = axons(curQuery.axonId);
    curAxon.nodesScaled = bsxfun(@times, ...
        curAxon.nodes(:, 1:3), p.raw.voxelSize);
    
    curQueries = queries( ...
        queries.axonId == curQuery.axonId ...
      & queries.chiasmaId == curQuery.chiasmaId, :);
    
   [~, curSummary] = ...
        connectEM.splitChiasmataMulti( ...
            p, curAxon, curQueries, 'dryRun', true);
        
    curFound = curSummary.tracings{1};
    curFound = curFound.overlaps{curQuery.exitId};
    curFound = setdiff(curFound, curQuery.exitId);
    
    gtEval.id(curIdx) = curQuery.id;
    gtEval.found{curIdx} = setdiff(curFound, 0);
    gtEval.expected{curIdx} = setdiff(curQuery.gtExit, 0);
end

gtEval.correct = cellfun( ...
    @(e, f) isequal(e(:), f(:)), ...
    gtEval.expected, gtEval.found);

%%
fprintf('\n');
fprintf('Evaluation of attachment:\n\n');
disp(gtEval);

fprintf('Correct: %.2f %%\n', 100 * mean(gtEval.correct));
fprintf('Too many attachments: %.2f %%\n', 100 * mean( ...
    cellfun(@numel, gtEval.found) > cellfun(@numel, gtEval.expected)));
fprintf('Too few attachments: %.2f %%\n', 100 * mean( ...
    cellfun(@numel, gtEval.found) < cellfun(@numel, gtEval.expected)));
%}

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
numEndings = cellfun(@(o) sum(o > 0), queries.overlaps);

flightEval = table;
flightEval.twoEndings = (numEndings == 2);
flightEval.tooFewEndings = (numEndings < 2);

temp = varfun(@sum, flightEval);
temp.Properties.VariableNames = flightEval.Properties.VariableNames;

fprintf('\n');
fprintf('Flight path evaluation:\n\n');
disp(temp);

%% Look at flight paths with too few attachments
%{
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
%}

%% calculate partitioning
[~, ~, queries.uniChiasmaId] = unique( ...
    queries(:, {'axonId', 'chiasmaId'}), 'rows');
chiasmaCount = max(queries.uniChiasmaId);

chiasma = table;
chiasma.lut = cell(chiasmaCount, 1);
chiasma.partition = cell(chiasmaCount, 1);
chiasma.valid = false(chiasmaCount, 1);

for curIdx = 1:chiasmaCount
    curQueries = queries(queries.uniChiasmaId == curIdx, :);
    curNrExits = size(curQueries, 1);
    
    curEdgesAll = curQueries.overlaps;
    curEdgesAll = cat(2, curEdgesAll{:})';
     
    curEdges = sort(curEdgesAll, 2);
    curEdges(~all(curEdges, 2), :) = [];
    curEdges = unique(curEdges, 'rows');
    curEdges = reshape(curEdges, [], 2);
    
    % Find parition
    curAdj = sparse( ...
        curEdges(:, 2), curEdges(:, 1), ...
        true, curNrExits, curNrExits);
    
   [~, curLUT] = graphconncomp(curAdj, 'Directed', false);
    curLUT = reshape(curLUT, [], 1);
   
    % Build partition
    curPartition = accumarray(curLUT, 1);
    curPartition = sort(curPartition, 'descend');
    
    % For each component with at least two elements, there exists an edge.
    % Hence, there is no excuse of the flight paths involved in these
    % components to be dangling.
    curIsValid = ~any( ...
        curEdgesAll(:, 2) == 0 & ismember( ...
        curEdgesAll(:, 1), curEdgesAll(:, 2)));
    
    chiasma.lut{curIdx} = curLUT;
    chiasma.partition{curIdx} = curPartition;
    chiasma.valid(curIdx) = curIsValid;
end

%% evaluate partitions
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

%% look at partitions
%{
thisValid = true;
thisPartition = [3; 1];

thisPartitionStr = {'invalid', 'valid'};
thisPartitionStr = strcat( ...
    thisPartitionStr{1 + thisValid}, ...
    '-', makePartitionStr(thisPartition));

chiasmaIds = find(arrayfun( ...
    @(v, p) v == thisValid && ...
    isequal(p{1}, thisPartition), ...
    chiasma.valid, chiasma.partition));
chiasmaIds = reshape(chiasmaIds, 1, []);

rng(0);
randIds = randperm(numel(chiasmaIds));
randIds = chiasmaIds(randIds(1:min(10, numel(chiasmaIds))));

curOutputDir = fullfile( ...
    outputDir, sprintf('%s-partition', thisPartitionStr));
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
        '%s-partition-%d.nml', thisPartitionStr, curIdx)));
end
%}

%% select chiasmata to split and flight paths to apply
thisValid = true;
thisPartition = [2; 2];

chiasmaIds = find(arrayfun( ...
    @(v, p) v == thisValid && ...
    isequal(p{1}, thisPartition), ...
    chiasma.valid, chiasma.partition));
chiasmaIds = reshape(chiasmaIds, 1, []);

theseQueries = queries( ...
    ismember(queries.uniChiasmaId, chiasmaIds), :);
theseQueries.execute = false(size(theseQueries, 1), 1);

% find minimal set of flight paths to execute
for curIdx = chiasmaIds
    curQueryIds = find(theseQueries.uniChiasmaId == curIdx);
    curQueries = theseQueries(curQueryIds, :);
    
    curNrExits = size(curQueries, 1);
    curQueryEdges = cell2mat(curQueries.overlaps(:)');
    curQueryEdges = sort(curQueryEdges, 1)';
    
    curAdj = curQueryEdges;
    curAdj(~all(curAdj, 2), :) = [];
    
    % build minimal spanning tree per component
    curAdj = sparse(curAdj(:, 2), curAdj(:, 1), 1, curNrExits, curNrExits);
   	curAdj = graphminspantree(curAdj, 'Method', 'Kruskal');
    
    clear curEdges;
   [curEdges(:, 2), curEdges(:, 1)] = find(curAdj);
   
    % find flights which form the minimal spanning tree
   [~, curExecIds] = ismember(curEdges, curQueryEdges, 'rows');
    theseQueries.execute(curQueryIds(curExecIds)) = true;
end