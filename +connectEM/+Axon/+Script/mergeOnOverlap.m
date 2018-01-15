% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');

% for prototyping
% aggloIds = [2121, 12348]; % example #1 works
% aggloIds = [2252, 14719]; % example #2 works
% aggloIds = [ 303,  5567]; % example #3 works

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile);
axons = axons.axons(axons.indBigAxons);

%%
voxelSize = param.raw.voxelSize;
nodesA = axons(aggloIds(1)).nodes(:, 1:3) .* voxelSize;
nodesB = axons(aggloIds(2)).nodes(:, 1:3) .* voxelSize;

lenA = pathLen(nodesA);
lenB = pathLen(nodesB);

% make sure that `lenA > lenB`
if lenA < lenB
    aggloIds = fliplr(aggloIds);
    
    temp = nodesA;
    nodesA = nodesB;
    nodesB = temp;
    clear temp;
    
    temp = lenA;
    lenA = lenB;
    lenB = temp;
    clear temp;
end

% calculate pair-wise distance
distMat = pdist2(nodesA, nodesB);
distVec = reshape(min(distMat, [], 1), [], 1);

% discard common nodes
nodeIdsB = find(distVec > 100);
nodesB = nodesB(nodeIdsB, :);
distMat = distMat(:, nodeIdsB);
distVec = distVec(nodeIdsB);

% discard edges with discarded nodes
edgesB = axons(aggloIds(2)).edges;
[~, edgesB] = ismember(edgesB, nodeIdsB);
edgesB = edgesB(all(edgesB, 2), :);
assert(issorted(edgesB, 2));

% calculate connected components
adjMatB = sparse( ...
    edgesB(:, 2), edgesB(:, 1), 1, ...
    size(nodesB, 1), size(nodesB, 1));
[numCompsB, lutB] = graphconncomp( ...
    adjMatB, 'Directed', false);

edgesAB = nan(numCompsB, 2);
for curIdx = 1:numCompsB
    curNodeIds = find(lutB == curIdx);
    curLen = pathLen(nodesB(curNodeIds, :));
    
    % ignore tiny components
    if curLen < 2000; continue; end
    
    curDistMat = distMat(:, curNodeIds);
   [~, curMinIdx] = min(curDistMat(:));
   [curNodeA, curNodeB] = ind2sub( ...
       size(curDistMat), curMinIdx);
   
    edgesAB(curIdx, 1) = curNodeA;
    edgesAB(curIdx, 2) = curNodeIds(curNodeB);
end

compIdsB = find(~isnan(edgesAB(:, 1)));
nodeIdsB = find(ismember(lutB, compIdsB));
nodesB = nodesB(nodeIdsB, :);

edgesAB = edgesAB(compIdsB, :);
[~, edgesAB(:, 2)] = ismember(edgesAB(:, 2), nodeIdsB);
[~, edgesB] = ismember(edgesB, nodeIdsB);
edgesB = edgesB(all(edgesB, 2), :);

% globalize edges
edgesB = edgesB + size(nodesA, 1);
edgesAB(:, 2) = edgesAB(:, 2) + size(nodesA, 1);

nodesAB = cat(1, nodesA, nodesB);
edgesA = axons(aggloIds(1)).edges;
edgesAB = cat(1, edgesA, edgesAB, edgesB);

skel = skeleton();
skel = skel.addTree( ...
    'Original A', ...
    axons(aggloIds(1)).nodes(:, 1:3), ...
    axons(aggloIds(1)).edges);
skel = skel.addTree( ...
    'Original B', ...
    axons(aggloIds(2)).nodes(:, 1:3), ...
    axons(aggloIds(2)).edges);
skel = skel.addTree( ...
    'Merged A+B', ...
    round(nodesAB ./ voxelSize), ...
    edgesAB);

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/merged.nml');

function len = pathLen(nodesNm)
    if size(nodesNm, 1) < 2
        len = 0;
        return;
    end
    
    adj = squareform(pdist(nodesNm));
    mst = graphminspantree(sparse(adj), 'Method', 'Kruskal');
    len = full(sum(mst(:)));
end