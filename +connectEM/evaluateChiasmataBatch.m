% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

clear;
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

fprintf('Flight path evaluation:\n\n');
disp(temp);

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