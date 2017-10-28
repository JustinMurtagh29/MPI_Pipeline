% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_08_a.mat');
outputDir = '/home/amotta/Desktop/invalid-tasks';

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

data = load(axonFile, 'info', 'oldAxons');
oldAxons = data.oldAxons.axons;
dataFiles = data.info.param.dataFiles;
clear data;

% load flights used in axons
dataFileCount = numel(dataFiles);
queries = cell(dataFileCount, 1);
flights = cell(dataFileCount, 1);

for curFileIdx = 1:numel(dataFiles)
    curDataFile = dataFiles{curFileIdx};
    curData = load(curDataFile);
    
    % original queries
    curQueries = table;
    curQueries.axonId = curData.queries(:, 1);
    curQueries.chiasmaId = curData.queries(:, 2);
    curQueries.exitId = curData.queries(:, 3);
    curQueries.exitNodeId = curData.queries(:, 8);
    curQueries.seedPos = curData.queries(:, 4:6);
    curQueries.centerNodeId = curData.queries(:, 7);
    curQueries.taskId = curData.taskIds;
    curFlights = curData.ff;
    clear curData;
    
    queries{curFileIdx} = curQueries;
    flights{curFileIdx} = curFlights;
    clear curQueries curFlights;
end

queries = cat(1, queries{:});
flights = Util.concatStructs('last', flights{:});

%% restrict to flights with comment
hasComment = ~cellfun(@isempty, flights.comments);
fprintf('# queries: %d\n', numel(hasComment));
fprintf('# queries with comment: %d\n', sum(hasComment));

flights = structfun( ...
    @(vals) vals(hasComment), ...
    flights, 'UniformOutput', false);

[~, queries.flightRow] = ismember( ...
    queries.taskId, flights.filenamesShort);
queries(~queries.flightRow, :) = [];

%% show random examples
rng(0);
rows = randperm(size(queries, 1));

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for curIdx = 4
    curRow = rows(curIdx);
    curQuery = queries(curRow, :);
    curAxon = oldAxons(curQuery.axonId);
    
    curComments = cell(size(curAxon.nodes, 1), 1);
    curComments(:) = {''};
    curComments{curQuery.centerNodeId} = 'Chiasma';
    curComments{curQuery.exitNodeId} = 'Query seed';
    
    curSkel = skeleton();
    curSkel = curSkel.addTree( ...
        'Axon', curAxon.nodes(:, 1:3), ...
        curAxon.edges, [], [], curComments);
    
    curFlightNodes = flights.nodes{curQuery.flightRow};
    curSkel = curSkel.addTree('Flight', curFlightNodes);
    
    curSkelName = sprintf( ...
        '%03d_axon-%d_node-%d.nml', ...
        curIdx, curQuery.axonId, curQuery.centerNodeId);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel.write(fullfile(outputDir, curSkelName));
end