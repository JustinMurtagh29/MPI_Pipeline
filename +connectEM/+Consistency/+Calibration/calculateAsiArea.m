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
assert(numel(trees.z) == numel(unique(trees.z)));

maxNodeCount = max(cellfun(@(n) numel(n.id), trees.nodes));

%% Generate vertices for rough mesh
% NOTE(amotta): `curT` and `curB` are the top and bottom skeletons
clear cur*;
curT = struct;

mesh = struct;
mesh.nodes = cell(height(trees), 1);

for curIdB = 1:height(trees)
    curB = struct;
    curB.nodes = trees.nodes{curIdB};
    curB.edges = trees.edges{curIdB};
    
    % Convert to tree-local node indices
    curB.edges = [curB.edges.source, curB.edges.target];
   [~, curB.edges] = ismember(curB.edges, curB.nodes.id);
    curB.edges = sort(curB.edges, 2, 'ascend');
    
    curB.nodes = [curB.nodes.x, curB.nodes.y, curB.nodes.z];
    
    % Make sure that path is linear
    curNodeOccur = accumarray(curB.edges(:), 1, [size(curB.nodes, 1), 1]);
    curEndIds = find(curNodeOccur == 1);
    assert(isempty(setdiff(curNodeOccur, 1:2)));
    assert(numel(curEndIds) == 2);
    
    if curIdB > 1
        % NOTE(amotta): Each linear path has two ends, d'uh! In the first
        % Z-slice we pick an arbitrary one of the two endings as the start
        % of the linear path. But starting from the second Z-slice we want
        % to pick the ending which is closer to the start of the path in
        % the previous Z-slice.
        curDist = pdist2( ...
            curT.nodes(1, :) .* voxelSize, ...
            curB.nodes(curEndIds, :) .* voxelSize);
        
       [~, curSortIds] = sort(curDist, 'ascend');
        curEndIds = curEndIds(curSortIds);
    end
    
    % Sort nodes in order of linear path
    curSortIds = graph(curB.edges(:, 1), curB.edges(:, 2));
    curSortIds = curSortIds.shortestpath(curEndIds(1), curEndIds(2));
    
    curB.nodes = curB.nodes(curSortIds, :);
    curB = rmfield(curB, 'edges');
    
    % NOTE(amotta): Interpolate linear path. We're not doing this to smooth
    % the curves or get better path length. It's just soooo much easier to
    % generate a quad mesh if both tracings have the same number of nodes.
    curLens = curB.nodes(2:end, :) - curB.nodes(1:(end - 1), :);
    curLens = [-eps; cumsum(sqrt(sum((curLens .* voxelSize) .^ 2, 2)))];
    curLens = curLens / curLens(end);
    
    curAlpha = reshape(linspace(0, 1, maxNodeCount), [], 1);
    curNodeIds = arrayfun(@(p) find(curLens < p, 1, 'last'), curAlpha);
    
    curAlpha = ...
        (curAlpha - curLens(curNodeIds)) ...
     ./ (curLens(curNodeIds + 1) - curLens(curNodeIds));
    curAlpha = 1 - curAlpha;
    
    curB.nodes = ...
        curAlpha .* curB.nodes(curNodeIds, :) ...
      + (1 - curAlpha) .* curB.nodes(curNodeIds + 1, :);
  
    % Save nodes and move on
    mesh.nodes{curIdB} = curB.nodes;
    curT = curB;
end

mesh.nodes = cell2mat(mesh.nodes);

%% Build edges of triangle mesh
clear cur*;

% Edges within linear path
curIds = reshape(1:(maxNodeCount - 1), [], 1);
curIds = curIds + ((1:height(trees)) - 1) * maxNodeCount;
curEdgesXY = [curIds(:), curIds(:) + 1];

% Edges to partner node in linear path below
curIds = reshape(1:maxNodeCount, [], 1);
curIds = curIds + ((1:(height(trees) - 1)) - 1) * maxNodeCount;
curEdgesZ = [curIds(:), curIds(:) + maxNodeCount];

% Diagonal edges to next node after partner node in linear path below
curIds = reshape(1:(maxNodeCount - 1), [], 1);
curIds = curIds + ((1:(height(trees) - 1)) - 1) * maxNodeCount;
curEdgesXYZ = [curIds(:), curIds(:) + maxNodeCount];

mesh.edges = [curEdgesXY; curEdgesZ; curEdgesXYZ];

%% Smooth mesh
clear cur*;

% NOTE(amotta); Build system of equations
% * Part A of the coefficient matrix represents the edges of the mesh. We
%   try to smooth the surface by minimizing the squared distance to all
%   neighboring nodes.
% * Part B of the coefficient matrix represents the manually placed nodes.
%   In the absence of part B, all nodes would collapse to a single point
%   during least-squares optimization. Instead, we add "virtual" edges to
%   the original nodes.
curNodeIds = [ ...
    (1:size(mesh.edges, 1))', mesh.edges(:, 1); ...
    (1:size(mesh.edges, 1))', mesh.edges(:, 2)];
curEntries = [ ...
    +ones(size(mesh.edges, 1), 1); ...
    -ones(size(mesh.edges, 1), 1)];
curMatA = sparse( ...
    curNodeIds(:, 1), curNodeIds(:, 2), curEntries, ...
    size(mesh.edges, 1), size(mesh.nodes, 1));

curNodeIds = [ ...
    (1:size(mesh.nodes, 1))', (1:size(mesh.nodes, 1))'];
curEntries = [ ...
    +ones(size(mesh.nodes, 1), 1)];
curMatB = sparse( ...
    curNodeIds(:, 1), curNodeIds(:, 2), curEntries, ...
    size(mesh.nodes, 1), size(mesh.nodes, 1));

% NOTE(amotta): All edge constraints (part A of coefficient matrix) have
% unit weight. To preserve the circumference of the axon-spine interface
% while allowing more smoothing in its center, we weight the virtual edges
% inversely proportional to the number of "actual edges". Nodes at the
% circumference have few "actual neighbors", which makes sure that the
% virtual edge weight is high, keeping them close to home.
curWeights = [ ...
    ones(size(mesh.edges, 1), 1); ...
    2 * 6 ./ accumarray(mesh.edges(:), 1)];
curMat = [curMatA; curMatB] .* curWeights;

mesh.nodesSmooth = mesh.nodes;
for curDimId = 1:size(mesh.nodes, 2)
    curTarget = [zeros(size(mesh.edges, 1), 1); mesh.nodes(:, curDimId)];
    curNodesSmooth = full(curMat \ (curTarget .* curWeights));
    mesh.nodesSmooth(:, curDimId) = curNodesSmooth;
end

%% Debugging: Show isosurface
clear cur*;

% Generate triangles
curFaces = reshape(1:(maxNodeCount - 1), [], 1);
curFaces = curFaces + ((1:(height(trees) - 1)) - 1) * maxNodeCount;
curFaces = reshape(curFaces, [], 1);

curFaces = [ ...
    (curFaces + [0, 1, (maxNodeCount + 1)]);
    (curFaces + [0, maxNodeCount + [0, 1]])];
curFaces = sort(curFaces, 2);

curFig = figure();
curAx = axes(curFig);
axis(curAx, 'equal');
view(curAx, 3);

origTriMesh = struct('vertices', mesh.nodes, 'faces', curFaces);
curOrigPatch = patch(origTriMesh);
curOrigPatch.FaceColor = 'blue';
curOrigPatch.EdgeColor = 'none';
curOrigPatch.FaceAlpha = 0.1;

smoothTriMesh = struct('vertices', mesh.nodesSmooth, 'faces', curFaces);
curSmoothPatch = patch(smoothTriMesh);
curSmoothPatch.FaceColor = 'red';
curOrigPatch.FaceAlpha = 0.5;

camlight(curAx);

%% Calculate areas
clear cur*;

origArea = origTriMesh.vertices .* voxelSize;
origArea = cross( ...
    origArea(origTriMesh.faces(:, 2), :) ...
  - origArea(origTriMesh.faces(:, 1), :), ...
    origArea(origTriMesh.faces(:, 3), :) ...
  - origArea(origTriMesh.faces(:, 1), :), 2);
origArea = sum(sqrt(sum(origArea .* origArea, 2)) / 2) / 1E6 %#ok

smoothArea = smoothTriMesh.vertices .* voxelSize;
smoothArea = cross( ...
    smoothArea(smoothTriMesh.faces(:, 2), :) ...
  - smoothArea(smoothTriMesh.faces(:, 1), :), ...
    smoothArea(smoothTriMesh.faces(:, 3), :) ...
  - smoothArea(smoothTriMesh.faces(:, 1), :), 2);
smoothArea = sum(sqrt(sum(smoothArea .* smoothArea, 2)) / 2) / 1E6 %#ok
