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
clear mstAdj;

% find list of branch point candidates
nodeDegree = accumarray(mstEdges(:), 1, size(segIds));
tipNodeIds = find(nodeDegree == 1);
bpNodeIds = find(nodeDegree > 2);

%% alternative approach
[~, maxIdx] = max(reshape(mstDists(tipNodeIds, tipNodeIds), [], 1));
[~, maxIdx] = ind2sub([numel(tipNodeIds), numel(tipNodeIds)], maxIdx);

refId = tipNodeIds(maxIdx);
otherIds = setdiff(tipNodeIds, refId);

% orient edges towards reference node
flipMask = ...
    mstDists(mstEdges(:, 2), refId) ...
  > mstDists(mstEdges(:, 1), refId);
mstEdges(flipMask, :) = ...
    fliplr(mstEdges(flipMask, :));
% mstEdges(:, 1) is source node
% mstEdges(:, 2) is destination node

% determine predecessor for each node
nodePredLUT = zeros(numel(segIds), 3);
nodePredLUT(mstEdges(:, 1), 1) = mstEdges(:, 2);
nodePredLUT(mstEdges(:, 1), 2) = 1:size(mstEdges, 1);
nodePredLUT(mstEdges(:, 1), 3) = sqrt(sum((...
    segPoints(segIds(mstEdges(:, 1)), :) ...
  - segPoints(segIds(mstEdges(:, 2)), :)) .^ 2, 2));

edgeToTipDists = zeros(size(mstEdges, 1), 1);
for curSeedId = reshape(otherIds, 1, [])
    curNodeId = curSeedId;
    curTipDist = 0;
    
    while curNodeId ~= refId
        curNodeLUT = shiftdim(nodePredLUT(curNodeId, :));
        
        curTipDist = curTipDist + curNodeLUT(3);
        if curTipDist < edgeToTipDists(curNodeLUT(2)); break; end
        edgeToTipDists(curNodeLUT(2)) = curTipDist;
        curNodeId = curNodeLUT(1);
    end
end

nodeSuccLUT = accumarray( ...
    mstEdges(:, 2), 1:size(mstEdges, 1), ...
    size(segIds), @(rows) {reshape(rows, [], 1)});

% check each branch point
maxLabelId = 1;
nodeLabels = zeros(size(segIds));
nodeLabels(refId) = maxLabelId;
nodeStack = refId;
while ~isempty(nodeStack)
    curNodeId = nodeStack(end);
    nodeStack(end) = [];
    
    curLabelId = nodeLabels(curNodeId);
    curSuccEdgeIds = nodeSuccLUT{curNodeId};
    
    curSuccNodeIds = mstEdges(curSuccEdgeIds, 1);
    nodeStack = cat(1, nodeStack(:), curSuccNodeIds);
    
    curSuccDists = edgeToTipDists(curSuccEdgeIds);
    
   [~, curMaxId] = max(curSuccDists);
    curSuccDiffBranch = curSuccDists > 3000;
    curSuccDiffBranch(curMaxId) = false;
    
    curSuccLabels = zeros(size(curSuccNodeIds));
    curSuccLabels(:) = curLabelId;
    curSuccLabels(curSuccDiffBranch) = ...
        maxLabelId + (1:sum(curSuccDiffBranch));
    
    nodeLabels(curSuccNodeIds) = curSuccLabels;
    maxLabelId = maxLabelId + sum(curSuccDiffBranch);
end
assert(all(nodeLabels));

nodeSets = accumarray( ...
    reshape(nodeLabels, [], 1), ...
    reshape(1:numel(segIds), [], 1), ...
    [], @(ids) {ids(:)});

%% alternative approach
[~, maxIdx] = max(reshape(mstDists(tipNodeIds, tipNodeIds), [], 1));
[~, maxIdx] = ind2sub([numel(tipNodeIds), numel(tipNodeIds)], maxIdx);

refId = tipNodeIds(maxIdx);
otherIds = setdiff(tipNodeIds, refId);

% orient edges towards reference node
flipMask = ...
    mstDists(mstEdges(:, 2), refId) ...
  > mstDists(mstEdges(:, 1), refId);
mstEdges(flipMask, :) = ...
    fliplr(mstEdges(flipMask, :));

% determine predecessor for each node
nodePredLUT = zeros(size(segIds));
nodePredLUT(mstEdges(:, 1)) = mstEdges(:, 2);

% sort by distance to ref
distToRef = mstDists(otherIds, refId);
[distToRef, sortIds] = sort(distToRef, 'descend');
otherIds = otherIds(sortIds);

expPairDists = squareform(pdist(distToRef));
actPairDists = mstDists(otherIds, otherIds);
delta = actPairDists - expPairDists;

otherGroups = zeros(size(otherIds));
while true
    curGroupId = find(~otherGroups, 1);
    if isempty(curGroupId); break; end
    
    % find other ids on same branch
    curMask = delta(:, curGroupId) < 7500;
    
    curGroupId = max(otherGroups) + 1;
    otherGroups(curMask) = curGroupId;
end

allGroups = zeros(size(segIds));
keyboard
for curGroupId = 1:max(otherGroups)
    curSeedIds = otherIds(otherGroups == curGroupId);
    curSeedIds = reshape(curSeedIds, 1, []);
    
    for curSeedId = curSeedIds
        curNodeId = curSeedId;
        while curNodeId > 0 && ...
                ~allGroups(curNodeId)
            allGroups(curNodeId) = curGroupId;
            curNodeId = nodePredLUT(curNodeId);
        end
    end
end
assert(all(allGroups));

nodeSets = accumarray( ...
    reshape(allGroups, [], 1), ...
    reshape(1:numel(segIds), [], 1), ...
    [], @(ids) {unique(cat(1, refId, ids(:)))});

%%
nodeSets = {reshape(1:numel(segIds), [], 1)};

curSetIdx = 1;
while curSetIdx <= numel(nodeSets)
    curSetNodeIds = nodeSets{curSetIdx};
    curSetSplit = false;
    
    for curNodeIdx = 1:numel(curSetNodeIds)
        curNodeId = curSetNodeIds(curNodeIdx);
        if ~bpCandMask(curNodeId); continue; end
        
        curNodeIds = pairDists(curSetNodeIds, curNodeId);
        curNodeIds = curSetNodeIds(curNodeIds < 4500);
        if isscalar(curNodeIds); continue; end
        
        % calculate expected distance
        curDiff = squareform(pdist(mstDists(curNodeIds, curNodeId)));
        curDiff = mstDists(curNodeIds, curNodeIds) - curDiff;
        
        % check if we neeed to split
       [curMaxDiff, curMaxIdx] = max(curDiff(:));
        % if curNodeId == 419; keyboard; end
        
        % pos #1: 8179 nm (for axon 10321, node 295)
        % pos #2: 8463 nm (for axon 10321, node 597)
        % pos #3: 8471 nm (for axon 10321, node 205)
        % neg #1: 6812 nm (for axon 10321, node 384) 
        if curMaxDiff < 7500; continue; end
        
       [curMaxRow, curMaxCol] = ...
            ind2sub(size(curDiff), curMaxIdx);
        curSplitNodeIds = curNodeIds([curMaxCol, curMaxRow]);
        
        % nodes need to be split into two branches
        curOtherIds = setdiff(curSetNodeIds, curNodeId);
        curOtherDists = mstDists(curOtherIds, curSplitNodeIds);
        curMask = curOtherDists(:, 2) > curOtherDists(:, 1);
        
        nodeSets{curSetIdx} = cat(1, curNodeId, curOtherIds(~curMask));
        nodeSets{end + 1} = cat(1, curNodeId, curOtherIds(curMask)); %#ok
        
        curSetSplit = true;
        break;
    end
    
    if ~curSetSplit
        curSetIdx = curSetIdx + 1;
    end
end

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