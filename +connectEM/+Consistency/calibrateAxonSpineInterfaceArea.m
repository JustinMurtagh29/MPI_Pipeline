% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%%
voxelSize = [11.24, 11.24, 28];
nmlFile = '/mnt/mpibr/data/Personal/mottaa/L4/2019-02-08-ASI-Area-Calibration/asi-29969_axon-65943_spine-head-99310.nml';

nml = slurpNml(nmlFile);
trees = NML.buildTreeTable(nml);

% Make sure each tree is confined to a single Z slice
trees.z = cellfun(@(n) unique(n.z), trees.nodes, 'UniformOutput', false);
assert(all(cellfun(@isscalar, trees.z)));

% Sort trees by Z slice
trees.z = cell2mat(trees.z);
trees = sortrows(trees, 'z', 'ascend');

maxNodeCount = max(cellfun(@(n) numel(n.id), trees.nodes));

%%
% NOTE(amotta): T for "top" and "B" for bottom
clear cur*;
curT = [];

curOut = struct;
curOut.nodes = cell(height(trees), 1);

for curIdB = 1:height(trees)
    curB = struct;
    curB.nodes = trees.nodes{curIdB};
    
    curB.edges = trees.edges{curIdB};
    curB.edges = [curB.edges.source, curB.edges.target];
    [~, curB.edges] = ismember(curB.edges, curB.nodes.id);
    curB.edges = sort(curB.edges, 2, 'ascend');
    
    curB.nodes = [curB.nodes.x, curB.nodes.y, curB.nodes.z];
    
    % Make sure that path is linear. Sort nodes linearly
    curNodeOccur = accumarray(curB.edges(:), 1, [size(curB.nodes, 1), 1]);
    curEndIds = find(curNodeOccur == 1);
    assert(isempty(setdiff(curNodeOccur, 1:2)));
    assert(numel(curEndIds) == 2);
    
    if curIdB > 1
        % Check which end is closer to start in previous slice
        curDist = pdist2( ...
            curT.nodes(1, :) .* voxelSize, ...
            curB.nodes(curEndIds, :) .* voxelSize);
        
       [~, curSortIds] = sort(curDist, 'ascend');
        curEndIds = curEndIds(curSortIds);
    end
    
    curSortIds = graph(curB.edges(:, 1), curB.edges(:, 2));
    curSortIds = curSortIds.shortestpath(curEndIds(1), curEndIds(2));
    
    curB = rmfield(curB, 'edges');
    curB.nodes = curB.nodes(curSortIds, :);
    
    % Interpolate
    curLens = curB.nodes(2:end, :) - curB.nodes(1:(end - 1), :);
    curLens = [-eps; cumsum(sqrt(sum((curLens .* voxelSize) .^ 2, 2)))];
    curLens = curLens / curLens(end);
    
    curAlpha = reshape(linspace(0, 1, maxNodeCount), [], 1);
    curIds = arrayfun(@(p) find(curLens < p, 1, 'last'), curAlpha);
    curIds = [curIds, (curIds + 1)];
    
    curAlpha = ...
        (curAlpha - curLens(curIds(:, 1))) ...
     ./ (curLens(curIds(:, 2)) - curLens(curIds(:, 1)));
    curAlpha = 1 - curAlpha;
    
    curB.nodes = ...
        curAlpha .* curB.nodes(curIds(:, 1), :) ...
      + (1 - curAlpha) .* curB.nodes(curIds(:, 2), :);
  
    curOut.nodes{curIdB} = curB.nodes;
    
    curT = curB;
end

% Build triangles
curEdgesXY = (1:maxNodeCount * height(trees)) - 1;
curEdgesXY(1:maxNodeCount:end) = [];
curEdgesXY = reshape(curEdgesXY, [], 1);
curEdgesXY = [curEdgesXY, (curEdgesXY + 1)];

curEdgesZ = reshape(1:maxNodeCount, [], 1);
curEdgesZ = curEdgesZ + ((1:(height(trees) - 1)) - 1) * maxNodeCount;
curEdgesZ = curEdgesZ(:) + [0, maxNodeCount];

curEdgesXZ = reshape(1:(maxNodeCount - 1), [], 1);
curEdgesXZ = curEdgesXZ + ((1:(height(trees) - 1)) - 1) * maxNodeCount;
curEdgesXZ = curEdgesXZ(:) + [0, (maxNodeCount + 1)];

curOut.nodes = cell2mat(curOut.nodes);
curOut.edges = [curEdgesXY; curEdgesZ; curEdgesXZ];

% Smooth
curIds = [ ...
    (1:size(curOut.edges, 1))', curOut.edges(:, 1); ...
    (1:size(curOut.edges, 1))', curOut.edges(:, 2)];
curEntries = [ ...
    +ones(size(curOut.edges, 1), 1); ...
    -ones(size(curOut.edges, 1), 1)];
curMatA = sparse( ...
    curIds(:, 1), curIds(:, 2), curEntries, ...
    size(curOut.edges, 1), size(curOut.nodes, 1));

curIds = [ ...
    (1:size(curOut.nodes, 1))', (1:size(curOut.nodes, 1))'];
curEntries = [ ...
    +ones(size(curOut.nodes, 1), 1)];
curMatB = sparse( ...
    curIds(:, 1), curIds(:, 2), curEntries, ...
    size(curOut.nodes, 1), size(curOut.nodes, 1));

curWeights = [ ...
    ones(size(curOut.edges, 1), 1); ...
    9 ./ accumarray(curOut.edges(:), 1)];
curMat = [curMatA; curMatB] .* curWeights;

curOut.nodesNew = curOut.nodes;
for curDimId = 1:3
    curGoal = [zeros(size(curOut.edges, 1), 1); curOut.nodes(:, curDimId)];
    curOut.nodesNew(:, curDimId) = full(curMat \ (curGoal .* curWeights));
end

curFaces = reshape(1:(maxNodeCount - 1), [], 1);
curFaces = curFaces + ((1:(height(trees) - 1)) - 1) * maxNodeCount;
curFaces = reshape(curFaces, [], 1);

curFaces = [ ...
    (curFaces + [0, 1, (maxNodeCount + 1)]);
    (curFaces + [0, maxNodeCount + [0, 1]])];

curP = struct;
curP.vertices = curOut.nodesNew;
curP.faces = sort(curFaces, 2);

figure;
axis equal
view(3);

curFigP = patch(curP);
curFigP.EdgeColor = 'none';
curFigP.FaceColor = 'red';
curFigP.FaceAlpha = 0.5;

curP.vertices = curOut.nodes;

curFigP = patch(curP);
curFigP.EdgeColor = 'none';
curFigP.FaceColor = 'blue';
curFigP.FaceAlpha = 0.5;

camlight
