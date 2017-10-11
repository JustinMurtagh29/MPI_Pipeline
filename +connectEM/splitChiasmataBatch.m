% This is an "improved" version of +connectEM/splitChiasmataMultiSuper
% which correctly handles the case where not all queries were answered yet.
%
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

%%
% Sort queries to make processing easier
queries = sortrows(queries, {'axonId', 'chiasmaId', 'exitId'});
flights = structfun(@(x) x(:), flights, 'UniformOutput', false);

% Assign flight paths to queries
[~, queries.flightId] = ismember( ...
    queries.taskId, flights.filenamesShort);

% Find chiasmata for which we have all queries answered
[~, ~, queries.uniChiasmaId] = ...
    unique(queries(:, {'axonId', 'chiasmaId'}), 'rows');
uniChiasmaDoneIds = find(accumarray( ...
    queries.uniChiasmaId, queries.flightId, [], @all));

% Limit ourselves to done chiasmata
queries = queries(ismember( ...
    queries.uniChiasmaId, uniChiasmaDoneIds), :);
queries.flightNodes = flights.nodes(queries.flightId);
queries.flightSegIds = flights.segIds(queries.flightId);

% Clean-up
queries(:, {'flightId', 'uniChiasmaId'}) = [];

%%
% Group queries by axons
[uniAxonIds, ~, queriesUniAxonIds] = unique(queries.axonId);
uniAxonQueries = accumarray( ...
    queriesUniAxonIds, (1:size(queries, 1))', ...
    [], @(rows) {queries(rows, :)});

% Start applying chiasmata
for curAxonIdx = 1:numel(uniAxonIds)
    curAxon = axons(uniAxonIds(curAxonIdx));
    curQueries = uniAxonQueries{curAxonIdx};
    
    connectEM.splitChiasmataMulti(p, curAxon, curQueries);
end