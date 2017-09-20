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

mergerAxon = axons(descIds(1));
otherAxons = axons(descIds(2:end));

%% inspect parts of the percolating merger
curBox = 5000;
curAxonId = descIds(1);

% restrict nodes and edges to bounding box
curAxon = mergerAxon;
curNodes = curAxon.nodes;
curEdges = curAxon.edges;

% make edges undirected
curEdges = sort(curEdges, 2);

% select random node
rng(3); % NOTE(amotta): increment seed from zero to 3
curNodeId = randperm(size(curNodes, 1), 1);
curNodeIdAbs = curNodeId;

% box centered on selected node
curBox = round(curBox ./ param.raw.voxelSize(:));
curBox = curBox * [-1, +1] / 2;
curBox = curBox + curNodes(curNodeId, 1:3)';

% reduce to bounding box
curIds = ...
    all(curNodes(:, 1:3) >= curBox(:, 1)', 2) ...
  & all(curNodes(:, 1:3) <= curBox(:, 2)', 2);
curNodes = curNodes(curIds, :);

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
[~, curEdges] = ismember(curEdges, curIds);
curEdges = curEdges(all(curEdges, 2), :);

% sanity check
assert(max(curEdges(:)) <= size(curNodes, 1));

%% show random (path-corrected) samples from other axons
curAxons = otherAxons;
curAxonLens = arrayfun( ...
    @(a) connectEM.getPathLengthFromNodes({a.nodes(:, 1:3)}), curAxons);

%%
rng(9); % NOTE(amotta): increment seed!
curAxonId = datasample(descIds(2:end), 1, 'Weights', curAxonLens);
curAxon = axons(curAxonId);

% fake entry
curNodeIdAbs = 0;

curNodes = curAxon.nodes;
curEdges = curAxon.edges;

%%
% NOTE(amotta): I should also split into connected components
% But let's see if webKNOSSOS will do this for us

curTreeName = sprintf('Axon %d, Node %d', curAxonId, curNodeIdAbs);
curFileName = sprintf('axon_%d__node_%d', curAxonId, curNodeIdAbs);
curFileName = fullfile('/home/amotta/Desktop', curFileName);

skel = skeleton();
skel = skel.addTree(curTreeName, curNodes(:, 1:3), curEdges);
skel = Skeleton.setParams4Pipeline(skel, param);
skel.write(curFileName);
