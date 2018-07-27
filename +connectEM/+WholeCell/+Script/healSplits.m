% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wholeCellFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

probThresh = 0.75;

info = Util.runInfo();
Util.showRunInfo(info);

%% Utility function
first = @(v) v(1);
wmean = @(v, w) sum(v .* (w / sum(w)));

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

graph = Graph.load(rootDir);
graph.id = reshape(1:height(graph), [], 1);

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

curMask = graph.borderIdx > 0;
graph.borderArea = nan(height(graph), 1);
graph.borderArea(curMask) = borderAreas(graph.borderIdx(curMask));
clear borderAreas;

maxSegId = Seg.Global.getMaxSegId(param);
[conn, syn] = connectEM.Connectome.load(param, connFile);

dendrites = load(conn.info.param.dendriteFile);
dendrites = dendrites.dendrites(conn.denMeta.parentId);

wholeCells = load(wholeCellFile);
wholeCells = wholeCells.wholeCells;

validDendAggloIds = find(~ismember( ...
    conn.denMeta.targetClass, {'Somata', 'WholeCell'}));
dendLUT = Agglo.buildLUT( ...
    maxSegId, conn.dendrites(validDendAggloIds), validDendAggloIds);

%%
clear cur*;
curCellId = 85;
curWholeCell = wholeCells(curCellId);
curSegIds = Agglo.fromSuperAgglo(curWholeCell);
curLUT = logical(Agglo.buildLUT(maxSegId, {curSegIds}));

curGraph = graph;
curGraph.cellMask = curLUT(curGraph.edges);
curGraph = curGraph(xor( ...
    curGraph.cellMask(:, 1), ...
    curGraph.cellMask(:, 2)), :);

curNeighDendIds = curGraph.edges(~curGraph.cellMask);
curNeighDendIds = setdiff(dendLUT(curNeighDendIds), 0);

curNeighT = table;
curNeighT.dendId = curNeighDendIds(:);
curNeighT.edgeIds = arrayfun( ...
    @(id) curGraph.id(any(ismember( ...
        curGraph.edges, conn.dendrites{id}), 2)), ...
    curNeighDendIds, 'UniformOutput', false);

curNeighT.meanEdgeProb = cellfun(@(ids) ...
    wmean(graph.prob(ids), graph.borderArea(ids)), curNeighT.edgeIds);

curNeighT = sortrows(curNeighT, 'meanEdgeProb', 'descend');
curNeighT(curNeighT.meanEdgeProb < probThresh, :) = [];

curNeighT.maxEdgeId = cellfun(@(ids) first( ...
    Util.sortBy(ids, graph.prob(ids), 'descend')), curNeighT.edgeIds);

curSkel = skeleton();
curSkel = Skeleton.setParams4Pipeline(curSkel, param);
curSkel = curSkel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

curSkel = Superagglos.toSkel(SuperAgglo.clean(curWholeCell), curSkel);
curSkel.names{end} = sprintf('Whole cell %d', curCellId);

curNumDigits = ceil(log10(1 + height(curNeighT)));
for curDendIdx = 1:height(curNeighT)
    curDendId = curNeighT.dendId(curDendIdx);
    curDendName = sprintf( ...
        '%0*d. Dendrite %d. Edge %d', ...
        curNumDigits, curDendIdx, curDendId, ...
        curNeighT.maxEdgeId(curDendIdx));
    
    curSkel = Superagglos.toSkel( ...
        SuperAgglo.clean(dendrites(curDendId)), curSkel);
    curSkel.names{end} = curDendName;
end
