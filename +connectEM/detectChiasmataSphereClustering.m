function output = detectChiasmata(p, nodesV, edges, visualize, outputFolder )
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here

% Create output folder if it does not exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Not very beautiful to temporary add to p structure
p.sphereRadiusOuter = 3500; % in nm
p.sphereRadiusInner = 3000; % in nm
p.minimumCosineDistance = 0.2;

% Scale to nm
nodes = bsxfun(@times, nodesV, p.raw.voxelSize);

% for each node with node degree > 2 ("marching sphere" approach to merger detection)
isIntersection = false(size(nodes,1),1);
nrExits = zeros(size(nodes,1),1);
for i=1:size(nodes,1)
    nodeDegree = sum(edges(:) == i);
    % If node degree of current node larger than 2
    if nodeDegree > 2
        % Step 1: Prune to sphere (around node in row i) and keep only CC of nodes connected to i
        [thisNodes, thisEdges, thisIdx, thisIdxOuter] = pruneToSphere(p, nodes, edges, i);
        % Sphere degree of node is defined as number of nodes connecting to the
        % outside of the outer sphere just pruned to
        edgesToOutside = thisEdges(any(ismember(thisEdges, thisIdxOuter),2),:);
        nodeSphereDegree = length(setdiff(edgesToOutside, thisIdxOuter));
        if nodeSphereDegree > 2
            % Step 2: Restrict to edges in outer sphere (and interpolate along them)
            thisNodesOuterSphere = sphereHull(p, thisNodes, thisEdges, thisIdx);
            % Step 3: Cluster into groups using single linkage clustering
            distances = squareform(pdist(thisNodesOuterSphere, 'cosine'));
            adjMatrix = sparse(double(distances < p.minimumCosineDistance));
            cc = Graph.findConnectedComponents(adjMatrix, false, false);
            if length(cc) > 2;
            	isIntersection(i) = true;
            end
            nrExits(i) = length(cc);
        else
            nrExits(i) = nodeSphereDegree;
        end
    end

    if visualize && isIntersection(i)
        figure('Position', [3841 1 1920 999]);
        % First subplot visualizing pruning to sphere (Step 1)
        subplot(2,2,1);
        visualizeSingleSphere(thisNodes, thisEdges, thisIdx, thisIdxOuter, p);
        title('Pruned to CC of edges within outer sphere');
        % Second subplot visualizing result of inital clustering
        subplot(2,2,2);
        visualizeClustering(thisNodes, thisEdges, thisIdx, thisIdxOuter, p, thisNodesOuterSphere, cc);
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
        title(['Final result: Detected intersection with ' num2str(nrExits(i)) ' exits']);
        pause(2);
        close all;
    end

end

% Find CC of detected intersections according to graph
if isempty(edges)
    cc = {};
else
    temp.edges = edges;
    cc = findCCaccordingToGraph(temp, find(isIntersection));
end
[~, centerOfCC] = cellfun(@(x)min(pdist2(bsxfun(@minus, nodes(x,:), mean(nodes(x,:),1)), [0 0 0])), cc);

% Find out where to query for each CC
queryIdx = cell(length(cc),1);
pos = cell(length(cc),1);
dir = cell(length(cc),1);
dist = cell(length(cc),1);
toDel = false(length(cc),1);
for i=1:length(cc)
    % Find all nodes outside CC of detected intersections that are neighbours to CC
    edgesIdx = any(ismember(edges, cc{i}),2);
    queryIdx{i} = setdiff(edges(edgesIdx,:), cc{i});
    % Keep edges not connected to CC
    edgesPruned = edges(~edgesIdx, :);
    nodeDegree = histc(edgesPruned(:), 1:size(nodes,1));
    % Find CC of graph with 5 or more elements after removing CC of intersections
    ccAfterPruning = Graph.findConnectedComponents(edgesPruned, false, true);
    ccAfterPruning = ccAfterPruning(cellfun('length', ccAfterPruning) > 4);
    % Keep only one one queryIdx for each of those CC (the one with maximum node degree)
    queryLocation = cellfun(@(x)intersect(x, queryIdx{i}), ccAfterPruning, 'uni', 0);
    queryLocation = queryLocation(~cellfun('isempty', queryLocation));
    [maxVal, maxIdx] = cellfun(@(x)max(nodeDegree(x)), queryLocation, 'uni', 0);
    if any(cellfun(@(x)x < 1, maxVal))
        % As Graph.findCC is used above with 3rd argument set to true
        error('This should not happen');
    end
    queryIdx{i} = cellfun(@(x,y)x(y), queryLocation, maxIdx);
    % Only generate queries for this CC if at least 2 queries
    if length(queryIdx{i}) > 1
        % Save position direction and length of query
        pos{i} = nodesV(queryIdx{i},:);
        dir{i} = bsxfun(@minus, pos{i}, nodesV(cc{i}(centerOfCC(i)),:));
        dist{i} = pdist2(bsxfun(@times, pos{i}, p.raw.voxelSize), bsxfun(@times, nodesV(cc{i}(centerOfCC(i)),:), p.raw.voxelSize), 'chebychev');
    else
        %toDel(i) = true;
    end
end
queryIdx(toDel) = [];
pos(toDel) = [];
dir(toDel) = [];
dist(toDel) = [];
cc(toDel) = [];
centerOfCC(toDel) = [];

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
output.distance = dist;

% Save result
save([outputFolder 'result.mat'], 'output');

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
