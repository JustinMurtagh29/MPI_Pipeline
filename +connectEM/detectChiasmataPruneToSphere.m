function [thisNodes, thisEdges, thisProb] = detectChiasmataPruneToSphere(nodes, edges, prob, p, i)
% Prune nodes, edges and prob to sphere around node in row i with param p

% Distance from all nodes to "current" node
assert(numel(i)==1);
thisDistance = pdist2(nodes(i,:), nodes);

% first make it small to keep it fast
thisNodeIdx = find(thisDistance < p.sphereRadiusOuter);
edges(~all(ismember(edges,thisNodeIdx),2),:) = [];
lookup(thisNodeIdx) = 1 : length(thisNodeIdx);
edges = lookup(edges);
nodes = nodes(thisNodeIdx, :);
i = lookup(i);

thisDistance = pdist2(nodes(i,:), nodes);
thisNodeIdx = thisDistance > p.sphereRadiusInner;
% rescue inner points that are not connected to center node within inner sphere
innerEdges = any(ismember(edges, find(~thisNodeIdx)),2);
innerConnected = Graph.findConnectedComponents(edges(innerEdges,:));
idx = find(cellfun(@(x)ismember(i,x),innerConnected));
for idx2 = setdiff(1: length(innerConnected), idx)
    thisNodeIdx(innerConnected{idx2}) = true;
end
%remove not connected skeletons
hereConnected = Graph.findConnectedComponents(edges);
idx = find(~cellfun(@(x)ismember(i,x),hereConnected));
idx = idx(:);
% assert(length(idx)<length(hereConnected));
for idx2 = idx'
    thisNodeIdx(hereConnected{idx2}) = false;
end


thisEdgeIdx = all(ismember(edges, find(thisNodeIdx)),2);
%this means that only nodes are kept that also have an edge outside of the sphere, should probably be changed
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
%visualizeSingleSphere(nodes, edges, prob, ~thisNodeIdx, thisNodeIdx, p);
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
