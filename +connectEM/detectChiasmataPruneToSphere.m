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
assert(length(idx)<length(hereConnected));
for idx2 = idx'
    thisNodeIdx(hereConnected{idx2}) = false;
end


% Keep all edges that have at least one node within outerSphere
thisEdgeIdx = all(ismember(edges, find(thisNodeIdx)),2);
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

end