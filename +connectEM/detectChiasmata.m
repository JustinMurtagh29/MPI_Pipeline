function output = detectChiasmata(nodesV, edges, visualize, outputFolder )
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here

% Create output folder if it does not exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Scale to nm
%nodes = bsxfun(@times, nodesV, [11.24 11.24 28]);
% Make sure edges are unique
[edges, ~, idxE] = unique(edges, 'rows');

% Some parameter for algorithm
p.sphereRadiusOuter = Inf; % in nm
p.sphereRadiusInner = 1000; % in nm
p.voxelSize = [11.24 11.24 28];
nodes=bsxfun(@times,nodesV,p.voxelSize);
% for each node ("marching sphere" approach to merger detection)
isIntersection = false(size(nodes,1),1);
nrExits = zeros(size(nodes,1),1);
for i=1:size(nodes,1)
    i
    [thisNodes, thisEdges, thisProb] = pruneToSphere(nodes, edges, ones(size(edges,1),1), p, i);
    conM = zeros(size(thisNodes, 1));
    conM(sub2ind(size(conM), thisEdges(:, 1), thisEdges(:, 2)))= 1;
    conM = conM+conM';
    [S,C] = graphconncomp(sparse(conM));
    if sum(arrayfun(@(idx)max(pdist2(thisNodes(C == idx, :), nodes(i,:))) > 4000, 1:S))>3
        isIntersection(i) = true;
        nrExits(i) = 4;
        % end
    end

    if false && isIntersection(i)
        figure('Position', [3841 1 1920 999]);
        % First subplot visualizing pruning to sphere (Step 1)
        subplot(2,2,1);
        visualizeSingleSphere(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter, p);
        colorbar;
        title('Pruned to CC of edges within outer sphere');
        % Second subplot visualizing result of inital clustering
        subplot(2,2,2);
        visualizeClustering(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter, p, thisNodes3, bestClustering.idxAll);
        colorbar;
        title('Clusters in sphere hull as detected by cluster visualization');
        % Third subplot
        subplot(2,2,3);
        imagesc(distances);
        axis equal; axis off;
        colorbar;
        title('Cosine distances between nodes in outer sphere');
        subplot(2,2,4);
        imagesc(distances > .1);
        axis equal; axis off;
        colorbar;
        title(['Final result: Detected intersection with ' num2str(nrExits(i)) ' exits']);
        % Save both as fig and png
        saveas(gcf, [outputFolder num2str(i) '.fig']);
        %img = getframe(gcf);
        %imwrite(img.cdata, [outputFolder num2str(i) '.png']);
        close all;
    end
end

% Find CC of detected intersections according to graph
temp.edges = edges;
cc = findCCaccordingToGraph(temp, find(isIntersection));
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
        dist{i} = pdist2(bsxfun(@times, pos{i}, p.voxelSize), bsxfun(@times, nodesV(cc{i}(centerOfCC(i)),:), p.voxelSize), 'chebychev');
    else
        toDel(i) = true;
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
output.prob = [];
output.isIntersection = isIntersection;
output.nrExits = nrExits;
output.ccNodeIdx = cc;
output.ccCenterIdx = cellfun(@(x,y)x(y), cc, mat2cell(centerOfCC, ones(size(centerOfCC,1),1), ones(size(centerOfCC,2),1)));
output.queryIdx = queryIdx;
output.position = pos;
output.direction = dir;
output.distance = dist;

if visualize
    % Write result to skletons for control (detection of intersections)
    comments = cell(size(nodesV,1), 1);
    comments(isIntersection) = arrayfun(@(x)strcat('intersection, ', num2str(x), ' exits'), nrExits(isIntersection), 'uni', 0);
    writeNml([outputFolder 'skelDetected.nml'], writeSkeletonFromNodesAndEdges({nodesV}, {edges}, {comments}, {['axon_' num2str(i, '%.5i')]}, {[0 0 1 1]}));
    % Write result to skletons for control (query visualization)
    clear comments;
    n{1} = nodesV(~isIntersection,:);
    % Keep only nodes & edges not part of any intersection
    edgeOffset = cumsum(accumarray(find(~isIntersection), 1))';
    e{1} = edges(~any(ismember(edges, find(isIntersection)),2),:);
    e{1} = edgeOffset(e{1});
    comments{1} = cell(size(n{1},1), 1);
    t{1} = 'original tree without detected intersections';
    col{1} = [0 0 1 1];
    for i=1:length(queryIdx)
        n{i+1} = cat(1, nodesV(output.ccCenterIdx(i),:), nodesV(output.queryIdx{i},:));
        e{i+1} = cat(2, ones(size(n{i+1},1)-1,1), (2:size(n{i+1},1))');
        comments{i+1} = cell(size(n{i+1},1), 1);
        comments{i+1}{1} = 'center node of intersection, edges inidcating queries';
        t{i+1} = ['intersection ' num2str(i)];
        col{i+1} = [1 0 0 1];
    end
    writeNml([outputFolder 'skelQueries.nml'], writeSkeletonFromNodesAndEdges(n, e, comments, t, col));
end

save([outputFolder 'result.mat'], 'output');

end

function [thisNodes, thisEdges, thisProb] = pruneToSphere(nodes, edges, prob, p, i)
% Prune nodes, edges and prob to sphere around node in row i with param p

% Distance from all nodes to "current" node
thisDistance = pdist2(nodes(i,:), nodes);

thisNodeIdx = thisDistance < p.sphereRadiusOuter & thisDistance > p.sphereRadiusInner;
%fix for not fully connected skeletons
hereConnected = Graph.findConnectedComponents(edges);
idx = cellfun(@(x)ismember(i,x),hereConnected);
if any(idx)
    thisNodeIdx = thisNodeIdx & ismember(1:length(thisDistance),hereConnected{idx})
else
    thisNodeIdx = false(size(thisNodeIdx));
end
% Keep all edges that have at least one node within outerSphere
thisEdgeIdx = any(ismember(edges, find(thisNodeIdx)),2);
thisIdxOuter = setdiff(unique(edges(thisEdgeIdx,:)), find(thisNodeIdx));
thisNodeIdx = unique(edges(thisEdgeIdx,:));
thisNodeIdx = thisNodeIdx(:);
thisOffset = cumsum(accumarray(thisNodeIdx, 1));
% Keep only the nodes, edges and prob within sphere
thisNodes = nodes(thisNodeIdx,:);
thisEdges = edges(thisEdgeIdx,:);
thisProb = prob(thisEdgeIdx);
% Rescale indices (thisEdges and "current" node (=i)) to be indices into
% thisNodes
thisEdges = thisOffset(thisEdges);
thisEdges = reshape(thisEdges,[],2); % Renumber according to new node indices
%thisIdx = thisOffset(i);
%thisIdxOuter = thisOffset(thisIdxOuter);

%[thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter] = restrictGraphToCC(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter);

end

function [thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter] = restrictGraphToCC(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter)

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
thisProb = thisProb(tempIdx);

end





function visualizeSingleSphere(nodes, edges, prob, currentNodeIdx, outerNodesIdx, p)

rP = p.sphereRadiusOuter;
rC = p.sphereRadiusInner;

% Define some values
nrBins = 10;
bins = linspace(0, 1, nrBins+1);

% RGB colors for indicating edge weight
edgeColors(:,1) = linspace(1,0,nrBins);
edgeColors(:,2) = linspace(0,1,nrBins);
edgeColors(:,3) = zeros(1,10);
% Calculate color for each edge
[row, col] = find(bsxfun(@gt, prob, bins(1:end-1)) & bsxfun(@le, prob, bins(2:end)));
[~, rIdx] = sort(row);
thisProbBinned = col(rIdx);
thisEdgeColors = edgeColors(thisProbBinned,:);

% Plot
hold on;
for i=1:size(edges,1)
    plot3([nodes(edges(i,1),1) nodes(edges(i,2),1)]', ...
        [nodes(edges(i,1),2) nodes(edges(i,2),2)]', ...
        [nodes(edges(i,1),3) nodes(edges(i,2),3)]', 'Color', thisEdgeColors(i,:));
end
plot3(nodes(currentNodeIdx, 1), nodes(currentNodeIdx, 2),nodes(currentNodeIdx, 3), 'xb', 'MarkerSize', 10);
plot3(nodes(outerNodesIdx, 1), nodes(outerNodesIdx, 2),nodes(outerNodesIdx, 3), 'xb', 'MarkerSize', 10);
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
colormap(edgeColors);
caxis([0 1]);

end

function visualizeClustering(nodes, edges, prob, currentNodeIdx, outerNodeIdx, p, clusterNodes, clusterIdx)

clusterNodes = bsxfun(@plus, clusterNodes, nodes(currentNodeIdx,:));

hold on;
visualizeSingleSphere(nodes, edges, prob, currentNodeIdx, outerNodeIdx, p);
clusters = unique(clusterIdx);
markers = {'or' 'og' 'ob' 'oy'};
for i=1:length(clusters)
    tempIdx = clusters(i) == clusterIdx;
    plot3(clusterNodes(tempIdx,1), clusterNodes(tempIdx,2), clusterNodes(tempIdx,3), markers{i}, 'MarkerSize', 10);
end

end
