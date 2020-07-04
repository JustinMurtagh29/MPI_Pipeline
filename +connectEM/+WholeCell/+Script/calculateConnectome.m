% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';

l4ConnRunId = '20190221T112510';
outDir = '/home/amotta/Desktop';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

curData = load(outputMapFile);
axonT = struct2table(curData.axonData, 'AsArray', true);

conn = curData.info.param.connFile;
[conn, syn] = connectEM.Connectome.load(param, conn);

[curDir, curFile] = fileparts(outputMapFile);
curAsiFile = sprintf('%s__%s_connectome.mat', curFile, l4ConnRunId);
curData = load(fullfile(curDir, curAsiFile));

synT = curData.synT;
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);

curAsiT = curData.asiT;
synT.type(:) = categorical({'Shaft'});
[~, curIds] = ismember(curAsiT.id, synT.id);
synT.type(curIds) = curAsiT.type;

%% Define utilities
clear cur*;

wmean = @(v, w) sum(v .* (w / sum(w)), 1);
segCom = @(ids) wmean(segPoints(ids, :), segSizes(ids));

%% Process axon NML files
axonT.somaPos = nan(height(axonT), 3);
axonT.cellId = nan(height(axonT), 1);

axonT.nodes(:) = {nan(0, 3)};
axonT.edges(:) = {nan(0, 2)};

for curId = 1:height(axonT)
    curNmlFile = axonT.nmlFile{curId};
    curSkel = skeleton(curNmlFile);
    
    curAxonTreeId = curSkel.getTreeWithName('Axon');
    curDendTreeId = curSkel.getTreeWithName('Dendrite');
    
    assert(isscalar(curAxonTreeId));
    assert(isscalar(curDendTreeId));
    
    curSomaNodeId = ...
        curSkel.getNodesWithComment( ...
            'Soma', curDendTreeId, 'insensitive');
    
	axonT.somaPos(curId, :) = ...
        curSkel.nodes{curDendTreeId}(curSomaNodeId, 1:3);
    
    axonT.nodes{curId} = curSkel.nodes{curAxonTreeId}(:, 1:3);
    axonT.edges{curId} = curSkel.edges{curAxonTreeId};
end

%% Find cell IDs for axon tracings
% NOTE(amotta): We will use these coordinates in conjuction with the soma
% nodes of the axon tracings to match tracings to cell IDs. This is so
% nasty that it hurts a bit...
clear cur*;

curMask = conn.denMeta.targetClass == 'Somata';
curSomaT = conn.denMeta(curMask, {'id', 'cellId'});
curSomaT = sortrows(curSomaT, 'cellId');

curSomaT.pos = conn.dendrites(curSomaT.id);
curSomaT.pos = cellfun(segCom, curSomaT.pos, 'UniformOutput', false);
curSomaT.pos = round(cell2mat(curSomaT.pos));

[~, curIds] = pdist2( ...
    curSomaT.pos .* param.raw.voxelSize, ...
    axonT.somaPos .* param.raw.voxelSize, ...
    'squaredeuclidean', 'Smallest', 1);

axonT.cellId = curSomaT.cellId(curIds);
assert(all(axonT.cellId));

assert(isequal( ...
    numel(axonT.cellId), ...
    numel(unique(axonT.cellId))));

[~, axonT.isInterneuron] = ismember(axonT.cellId, conn.denMeta.cellId);
axonT.isInterneuron = conn.denMeta.isInterneuron(axonT.isInterneuron);

axonT.numOutputSynapses = accumarray( ...
    synT.preAggloId, 1, [height(axonT), 1]);

%% Build connectivity matrix
% NOTE(amotta): Here, "cell" stands for soma-based neuron reconstructions.

cellSynT = synT;
cellSynT.postAggloId = conn.denMeta.cellId(cellSynT.postAggloId);
cellSynT = cellSynT(cellSynT.postAggloId > 0, :);

cellSynT.Properties.VariableNames{strcmp( ...
    cellSynT.Properties.VariableNames, 'preAggloId')} = 'preCellId';
cellSynT.Properties.VariableNames{strcmp( ...
    cellSynT.Properties.VariableNames, 'postAggloId')} = 'postCellId';

uniCellIds = union(cellSynT.preCellId, cellSynT.postCellId);
[~, curRelPreIds] = ismember(cellSynT.preCellId, uniCellIds);
[~, curRelPostIds] = ismember(cellSynT.postCellId, uniCellIds);

cellConn = accumarray( ...
    [curRelPreIds(:), curRelPostIds(:)], 1, ...
    [numel(uniCellIds), numel(uniCellIds)]);

%% Report
numAxons = height(axonT);
numInterneuronAxons = sum(axonT.isInterneuron);
totalAxonalPathLength = sum(axonT.pathLength) / 1E3;
averageAxonalPathLength = mean(axonT.pathLength) / 1E3;
medianAxonalPathLength = median(axonT.pathLength) / 1E3;

totalOutSynapses = sum(axonT.numOutputSynapses);
averageOutSynapses = mean(axonT.numOutputSynapses);
medianOutSynapses = median(axonT.numOutputSynapses);

numAxonsWithOutSynapses = sum(axonT.numOutputSynapses > 0);

numCellCellSynapses = sum(cellConn(:));
numCellCellConnections = sum(cellConn(:) > 0);
numCellsWithOutSynapsesOntoCells = numel(unique(cellSynT.preCellId));
numCellsWithInSynapsesFromCells = numel(unique(cellSynT.postCellId));

fprintf('Number of soma-based axons: %d\n', numAxons);
fprintf('  Number of interneuron axons: %d\n', numInterneuronAxons);
fprintf('Total axonal path length: %.3f mm\n', totalAxonalPathLength / 1E3);
fprintf('  Average: %.3f µm\n', averageAxonalPathLength);
fprintf('  Median: %.3f µm\n', medianAxonalPathLength);
fprintf('Total number of output synapses: %d\n', totalOutSynapses);
fprintf('  Average: %.2f\n', averageOutSynapses);
fprintf('  Median: %d\n', medianOutSynapses);
fprintf('Number of axons with output synapses: %d\n', numAxonsWithOutSynapses);
fprintf('\n');
fprintf('Number of cell-to-cell synapses: %d\n', numCellCellSynapses);
fprintf('Number of cell-to-cell connections: %d\n', numCellCellConnections);
fprintf('Number of cells with output synapses onto cells: %d\n', numCellsWithOutSynapsesOntoCells);
fprintf('Number of cells with intput synapses from cells: %d\n', numCellsWithInSynapsesFromCells);

% * Total path length
%   * median path length

