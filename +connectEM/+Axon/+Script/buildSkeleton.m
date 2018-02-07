% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

%% loading data
conn = load(connFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
segPoints = segPoints .* param.raw.voxelSize;

%% prototyping
axonId = 10321;
segIds = conn.axons{axonId};

%% prepare data
% build a skeleton representation
pairDists = squareform(pdist(segPoints(segIds, :)));

% calculate pair-wise distance (along skeleton)
mstAdj = graphminspantree(sparse(pairDists));
mstDists = graphallshortestpaths(mstAdj, 'Directed', false);

mstEdges = nan(nnz(mstAdj), 2);
[mstEdges(:, 2), mstEdges(:, 1)] = find(mstAdj);

% find list of branch point candidates
nodeDegree = accumarray(mstEdges(:), 1, size(segIds));
tipNodeIds = find(nodeDegree == 1);
bpNodeIds = find(nodeDegree > 2);

%% alternative approach
[~, maxIdx] = max(reshape(mstDists(tipNodeIds, tipNodeIds), [], 1));
[~, maxIdx] = ind2sub([numel(tipNodeIds), numel(tipNodeIds)], maxIdx);
refId = tipNodeIds(maxIdx);

[nodeLabels, nodeDists, mstEdges] = ...
    Graph.buildDendrogram(mstAdj, refId, 3500);

nodeSets = accumarray( ...
    reshape(nodeLabels, [], 1), ...
    reshape(1:numel(segIds), [], 1), ...
    [], @(ids) {ids(:)});

%% build skeleton
skelNodes = round(segPoints(segIds, :) ./ param.raw.voxelSize);
skelEdges = sort(mstEdges, 2);

skel = skeleton();
skel = skel.addTree(sprintf('Agglomerate %d', axonId), skelNodes, skelEdges);

for curSetIdx = 1:numel(nodeSets)
    curNodeIds = nodeSets{curSetIdx};
    curNodes = skelNodes(curNodeIds, :);
    
   [~, curEdges] = ismember(mstEdges, curNodeIds);
    curEdges = curEdges(all(curEdges, 2), :);
    curEdges = sort(curEdges, 2);
    
    curAdj = sparse( ...
        curEdges(:, 2), curEdges(:, 1), 1, ...
        numel(curNodeIds), numel(curNodeIds));
    assert( ...
        graphconncomp(curAdj, 'Directed', false) == 1, ...
        'Branch %d is not a connected component', curSetIdx);
    
    skel = skel.addTree( ...
        sprintf('Branch %d', curSetIdx), curNodes, curEdges);
end

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write(fullfile('/home/amotta/Desktop', sprintf('agglo-%d.nml', axonId)));

%% plot axonogram
branchEdges = nodeLabels(mstEdges);
branchEdges(~diff(branchEdges, 1, 2), :) = [];
branchEdges = unique(branchEdges, 'rows');

treeEdges = 1 + cat(1, [1, 0], branchEdges);

branchMaxDist = accumarray(nodeLabels, nodeDists, [], @max);
treeNodes = cat(2, branchMaxDist, zeros(size(branchMaxDist, 1), 2));
treeNodes = cat(1, zeros(1, 3), treeNodes);

tree = struct;
tree.name = sprintf('axon %d', axonId);
tree.edges = treeEdges;
tree.nodes = treeNodes;

axonogram(tree, [], 1, [1, 1, 1])