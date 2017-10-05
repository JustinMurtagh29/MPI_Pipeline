function output = detectChiasmataSphereClustering( ...
        p, nodesV, edges, visualize, outputFolder)
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here

if ~exist('outputFolder', 'var')
    outputFolder = [];
end

% Create output folder if it does not exist
if ~isempty(outputFolder) && ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Not very beautiful to temporary add to p structure
p.sphereRadiusOuter = 1500; % in nm
p.sphereRadiusInner = 1000; % in nm
p.minimumCosineDistance = 0.2;

% Scale to nm
nodes = bsxfun(@times, nodesV, p.raw.voxelSize);

% for each node with node degree > 2 ("marching sphere" approach to merger detection)
isIntersection = false(size(nodes,1),1);
nrExits = zeros(size(nodes,1),1);
for i=1:size(nodes,1)
    nrExits(i) = forNode(p, nodes, edges, i, visualize);
    isIntersection(i) = nrExits(i) > 3;
end

% Build chiasmata from chiasmatic nodes
[cc, centerOfCC] = ...
    connectEM.detectChiasmataNodesCluster(nodes, isIntersection);

% Find out where to query for each CC
queryIdx = cell(numel(cc), 1);
pos = cell(numel(cc),1);
dir = cell(numel(cc),1);
for i=1:numel(cc)
    [~, pos{i}, dir{i}, queryIdx{i}] = forNode( ...
        p, nodes, edges, cc{i}(centerOfCC(i)), false);
end

% Create an output structure
output.nodes = nodesV;
output.edges = edges;
output.isIntersection = isIntersection;
output.nrExits = nrExits;
output.ccNodeIdx = cc;
output.ccCenterIdx = cellfun(@(x,y)x(y), cc, num2cell(centerOfCC));
output.queryIdx = queryIdx;
output.position = pos;
output.direction = dir;

% Save result
if ~isempty(outputFolder)
    save(fullfile(outputFolder, 'result.mat'), 'output');
end

end

function [nrExits, pos, dir, queryIdx] = ...
        forNode(p, nodes, edges, i, visualize)
    nrExits = 0;
    pos = nan(0, 3);
    dir = nan(0, 3);
    queryIdx = nan(1, 0);
    
    %% analyse node vincinity
    nodeDegree = sum(edges(:) == i);
    if nodeDegree < 3; return; end;
    
    % Step 1: Prune to sphere (around node in row i) and keep only CC of nodes connected to i
    [thisNodes, thisEdges, thisIdx, thisIdxOuter] = pruneToSphere(p, nodes, edges, i);
    
    % Sphere degree of node is defined as number of nodes connecting to the
    % outside of the outer sphere just pruned to
    edgesToOutside = thisEdges(any(ismember(thisEdges, thisIdxOuter),2),:);
    nodeSphereDegree = length(setdiff(edgesToOutside, thisIdxOuter));
    if nodeSphereDegree < 3; return; end;
    
    % Step 2: Restrict to edges in outer sphere (and interpolate along them)
    hullNodes = sphereHull(p, thisNodes, thisEdges, thisIdx);
    
    % Step 3: Cluster into groups using single linkage clustering
    distances = squareform(pdist(hullNodes, 'cosine'));
    adjMatrix = sparse(double(distances < p.minimumCosineDistance));
    cc = Graph.findConnectedComponents(adjMatrix, false, false);
    
    exitIds = find(cellfun(@(idx) max(pdist2( ...
        hullNodes(idx, :), nodes(i, :))) > 3000, cc));
    cc = cc(exitIds);
    
    nrExits = numel(exitIds);
    if ~nrExits; return; end;
    
    %% determine queries, if desired
    if nargout > 2
        pos = nan(nrExits, 3);
        dir = nan(nrExits, 3);
        queryIdx = nan(1, nrExits);
        
        descale = @(nm) ceil(bsxfun(@times, nm, 1 ./ p.raw.voxelSize));
        
        for curIdx = 1:numel(cc)
            % find outer node corresponding to cluster
            curNodeIds = cc{curIdx};
            curNodePos = hullNodes(curNodeIds, :);
            
            % find node with largest distance from center
           [~, curMaxIdx] = max(sum(curNodePos .* curNodePos, 2));
            curNodePos = curNodePos(curMaxIdx, :);
            curNodePos = curNodePos + nodes(i, :);
            
            pos(curIdx, :) = descale(curNodePos);
            dir(curIdx, :) = pos(curIdx, :) - descale(nodes(i, :));
            
            % TODO(amotta): Implement!
            queryIdx(curIdx) = nan;
        end
    end
    
    %% possibly to visualization
    if visualize
        figure('Position', [3841 1 1920 999]);
        % First subplot visualizing pruning to sphere (Step 1)
        subplot(2,2,1);
        visualizeSingleSphere(thisNodes, thisEdges, thisIdx, thisIdxOuter, p);
        title('Pruned to CC of edges within outer sphere');
        % Second subplot visualizing result of inital clustering
        subplot(2,2,2);
        visualizeClustering(thisNodes, thisEdges, thisIdx, thisIdxOuter, p, hullNodes, cc);
        title('Clusters in sphere hull as detected by cluster visualization');
        % Third subplot
        subplot(2,2,3);
        imagesc(distances);
        axis equal; axis off;
        colorbar;
        title('Cosine distances between nodes in outer sphere');
        subplot(2,2,4);
        imagesc(distances > p.minimumCosineDistance);
        axis equal; axis off;
        colorbar;
        title(['Final result: Detected intersection with ' num2str(nrExits) ' exits']);
        pause(2);
        close all;
    end
end

function [thisNodes, thisEdges, thisIdx, thisIdxOuter] = pruneToSphere(p, nodes, edges, i)
% Prune nodes, edges to sphere of radius p.sphereRadiusOuter around node in row i
% This will keep all nodes that participate in an edge with one node closer
% than the given radius, returns new indices of current node (in former row
% i) and all indices to nodes that are outside of the sphere
% thisNodes will be centered on current node

% Distance from all nodes to "current" node
thisDistance = pdist2(nodes(i,:), nodes);
thisNodeIdx = thisDistance < p.sphereRadiusOuter;
% Keep all edges that have at least one node within outerSphere
thisEdgeIdx = any(ismember(edges, find(thisNodeIdx)),2);
thisIdxOuter = setdiff(unique(edges(thisEdgeIdx,:)), find(thisNodeIdx));
% Keep all nodes that are part of an edge in outerSphere
thisNodeIdx = unique(edges(thisEdgeIdx,:));
thisOffset = cumsum(accumarray(thisNodeIdx, 1));
% Keep only the nodes belong to an edge within sphere
thisNodes = nodes(thisNodeIdx,:);
thisEdges = edges(thisEdgeIdx,:);
% Rescale indices (thisEdges and "current" node (=i)) to be indices into
% thisNodes
thisEdges = thisOffset(thisEdges); % Renumber according to new node indices
thisIdx = thisOffset(i);
thisIdxOuter = thisOffset(thisIdxOuter);

% Center on current node
thisNodes = bsxfun(@minus, thisNodes, thisNodes(thisIdx,:));

[thisNodes, thisEdges, thisIdx, thisIdxOuter] = restrictGraphToCC(thisNodes, thisEdges, thisIdx, thisIdxOuter);

end

function [thisNodes, thisEdges, thisIdx, thisIdxOuter] = restrictGraphToCC(thisNodes, thisEdges, thisIdx, thisIdxOuter)

% Keep only CC connected to "current" node
thisCC = Graph.findConnectedComponents(thisEdges, false, false);
thisCC = thisCC{cellfun(@(x)any(ismember(x,thisIdx)), thisCC)};
thisNodes = thisNodes(thisCC,:);
thisIdx = find(ismember(thisCC,thisIdx));
thisIdxOuter = intersect(thisIdxOuter, thisCC);
for i=1:length(thisIdxOuter)
    thisIdxOuter(i) = find(ismember(thisCC,thisIdxOuter(i)));
end
tempIdx = all(ismember(thisEdges,thisCC),2);
thisEdges = thisEdges(tempIdx,:);
thisEdges = arrayfun(@(x)find(x == thisCC), thisEdges);

end


function thisNodes2 = sphereHull(p, nodes, edges, i)
% Remove all edges based on given probability threshold

% Distance of all nodes to current node
thisNodes = bsxfun(@minus, nodes, nodes(i,:));

% Interpolate
thisNodes2 = [];
for i=1:size(edges,1)
    startPos = thisNodes(edges(i,1),:);
    endPos = thisNodes(edges(i,2),:);
    lengthOfEdge = pdist2(startPos, endPos);
    % If necessary because there are duplicate nodes in ~50/40000 axons
    if lengthOfEdge
        nodesToAdd(:,1) = linspace(startPos(:,1), endPos(:,1), ceil(lengthOfEdge ./ 50));
        nodesToAdd(:,2) = linspace(startPos(:,2), endPos(:,2), ceil(lengthOfEdge ./ 50));
        nodesToAdd(:,3) = linspace(startPos(:,3), endPos(:,3), ceil(lengthOfEdge ./ 50));
        thisNodes2 = cat(1, thisNodes2, nodesToAdd);
    end
    clear nodesToAdd;
end

% Calculate spherical coordinates
[~, ~, r] = cart2sph(thisNodes2(:,1), thisNodes2(:,2), thisNodes2(:,3));

% Keep only nodes within given sphere hull
tempIdx = r > p.sphereRadiusInner & r < p.sphereRadiusOuter;
thisNodes2 = thisNodes2(tempIdx,:);

end


function [eva, best] = clusterOnUnitSphere(nodes)
% Cluster nodes using kmedoids with cosine distance, find optimal value
% of number clusters using silhouette criterion (for now)

% Define function used for calculation of distances between observations
options = statset('Display', 'off', 'MaxIter', 10);
fCluster = @(x,k)kmedoids(x, k, 'Algorithm', ...
    'pam', 'Distance', 'cosine', 'Options', options, ...
    'Replicates', 5, 'Start', 'sample');
eva = evalclusters(nodes, fCluster, 'silhouette', 'KList', 2:4, ...
    'Distance', 'cosine', 'ClusterPriors', 'equal');
% If more than one cluster is never detected
if ~isnan(eva.OptimalK)
    [best.idxAll, best.distCluster, ~, best.distNodes, best.idxMed] = fCluster(nodes, eva.OptimalK);
else
    best = [];
end

end

function visualizeSingleSphere(nodes, edges, currentNodeIdx, outerNodesIdx, p)

rP = p.sphereRadiusOuter;
rC = p.sphereRadiusInner;
edgeColor = 'g';

% Plot
hold on;
for i=1:size(edges,1)
    plot3([nodes(edges(i,1),1) nodes(edges(i,2),1)]', ...
        [nodes(edges(i,1),2) nodes(edges(i,2),2)]', ...
        [nodes(edges(i,1),3) nodes(edges(i,2),3)]', 'Color', edgeColor);
end
plot3(nodes(currentNodeIdx, 1), nodes(currentNodeIdx, 2),nodes(currentNodeIdx, 3), 'xb', 'MarkerSize', 10);
plot3(nodes(outerNodesIdx, 1), nodes(outerNodesIdx, 2),nodes(outerNodesIdx, 3), 'xr', 'MarkerSize', 10);
otherNodesIdx = setdiff(1:size(nodes,1), cat(1, currentNodeIdx, outerNodesIdx));
plot3(nodes(otherNodesIdx, 1), nodes(otherNodesIdx, 2),nodes(otherNodesIdx, 3), 'xy', 'MarkerSize', 10);
% Add spheres
[x,y,z] = sphere;
surf(rP*x+nodes(currentNodeIdx,1),rP*y+nodes(currentNodeIdx,2),rP*z+nodes(currentNodeIdx,3), 'EdgeColor', 'none', 'FaceColor', 'b');
surf(rC*x+nodes(currentNodeIdx,1),rC*y+nodes(currentNodeIdx,2),rC*z+nodes(currentNodeIdx,3), 'EdgeColor', 'none', 'FaceColor', 'b');
alpha(.2);
camlight;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

end

function visualizeClustering(nodes, edges, currentNodeIdx, outerNodeIdx, p, clusterNodes, cc)

clusterNodes = bsxfun(@plus, clusterNodes, nodes(currentNodeIdx,:));

hold on;
visualizeSingleSphere(nodes, edges, currentNodeIdx, outerNodeIdx, p);
markers = {'or' 'og' 'ob' 'oy'};
for i=1:length(cc)
    plot3(clusterNodes(cc{i},1), clusterNodes(cc{i},2), clusterNodes(cc{i},3), markers{i}, 'MarkerSize', 10);
end

end
