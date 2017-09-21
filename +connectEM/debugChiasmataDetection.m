% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% load parameters
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% load axons
rootDir = param.saveFolder;
data = load(fullfile(rootDir, 'aggloState', 'axons_05.mat'));

% only look at large enough axons
axons = data.axons(data.indBigAxons);

%% find percolating merger
% do this by counting segments
segmentCount = arrayfun(@(a) sum(a.nodes(:, end) > 0), axons);
[descSegmentCount, descIds] = sort(segmentCount, 'descend');

fprintf('Number of segments in largest agglomerates:\n');
disp(descSegmentCount(1:10));

%% inspect parts of the percolating merger
curBox = 5000;

curAxonId = descIds(1);
curAxon = axons(curAxonId);

% restrict nodes and edges to bounding box
curNodes = curAxon.nodes;
curEdges = curAxon.edges;

% keep track of node IDs
% TODO(amotta): use table
curNodeIds = 1:size(curNodes, 1);
curNodeIds = reshape(curNodeIds, [], 1);

% make edges undirected
curEdges = sort(curEdges, 2);

% select random node
rng(0); % NOTE(amotta): increment seed from zero to 3
curNodeId = randperm(size(curNodes, 1), 1);

% box centered on selected node
curBox = round(curBox ./ param.raw.voxelSize(:));
curBox = curBox * [-1, +1] / 2;
curBox = curBox + curNodes(curNodeId, 1:3)';

% reduce to bounding box
curIds = ...
    all(curNodes(:, 1:3) >= curBox(:, 1)', 2) ...
  & all(curNodes(:, 1:3) <= curBox(:, 2)', 2);
curNodes = curNodes(curIds, :);
curNodeIds = curNodeIds(curIds);

curIds = find(curIds);
[~, curNodeId] = ismember(curNodeId, curIds);
[~, curEdges] = ismember(curEdges, curIds);
curEdges = curEdges(all(curEdges, 2), :);

% sanity check
assert(max(curEdges(:)) <= size(curNodes, 1));

% build connected components
curNodeCount = size(curNodes, 1);
curAdjMat = sparse( ...
    curEdges(:, 2), curEdges(:, 1), 1, curNodeCount, curNodeCount);

[~, curCompIds] = graphconncomp(curAdjMat, 'Directed', false);
curCompId = curCompIds(curNodeId);

curIds = find(curCompIds == curCompId);
[~, curNodeId] = ismember(curNodeId, curIds);

curNodes = curNodes(curIds, :);
curNodeIds = curNodeIds(curIds);

[~, curEdges] = ismember(curEdges, curIds);
curEdges = curEdges(all(curEdges, 2), :);

% sanity check
assert(max(curEdges(:)) <= size(curNodes, 1));
assert(size(curNodeIds, 1) == size(curNodes, 1));

curNodeIdAbs = curNodeIds(curNodeId);

%% calculate path length for other axons
curAxonLens = arrayfun( ...
    @(a) connectEM.getPathLengthFromNodes({a.nodes(:, 1:3)}), ...
    axon(descIds(2:end)));

%% show random (path-corrected) samples from other axons
rng(9); % NOTE(amotta): increment seed!
curAxonId = datasample(descIds(2:end), 1, 'Weights', curAxonLens);
curAxon = axons(curAxonId);

curNodes = curAxon.nodes;
curEdges = curAxon.edges;
curNodeIdAbs = 0;

%% prepare to show chiasmatic nodes
% `connectEM.detectChiasmata` saves a `isIntersection` for each axon
% agglomerate. Let's use this to highlight nodes.

% load results of chiasmata detection
% see https://gitlab.mpcdf.mpg.de/connectomics/pipeline/blob/master/+connectEM/detectChiasmataSuper.m
chiasmDir = fullfile( ...
    param.saveFolder, 'chiasmata', ...
    sprintf('chiasmataX%d_%d', 32, floor(curAxonId / 100)), ...
    sprintf('visX%d_%d', 32, curAxonId));
chiasmFile = fullfile(chiasmDir, 'results.mat');

chiasmNodes = load(chiasmFile, 'output');
chiasmNodes = find(chiasmNodes.output.isIntersection);

% restrict to node subset
[~, chiasmNodes] = ismember(chiasmNodes, curNodeIds);
chiasmNodes = setdiff(chiasmNodes, 0);

%%
% NOTE(amotta): I should also split into connected components
% But let's see if webKNOSSOS will do this for us

curTreeName = sprintf('Axon %d, Node %d', curAxonId, curNodeIdAbs);
curFileName = sprintf('axon_%d__node_%d', curAxonId, curNodeIdAbs);
curFileName = fullfile('/home/amotta/Desktop', curFileName);

skel = skeleton();
skel = skel.addTree(curTreeName, curNodes(:, 1:3), curEdges);

% show chiasmata nodes in yellow
skel = skel.addNodesAsTrees(curNodes(chiasmaNodes, 1:3), 'Chiasma');
skel.colors((end - numel(chiasmaNodes) - 1):end) = {[1, 1, 0, 1]};

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write(curFileName);
