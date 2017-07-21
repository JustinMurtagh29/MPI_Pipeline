function [thisNodes, thisEdges, thisProb] = detectChiasmataPruneToSphere(nodes, edges, prob, p, i)
% Prune nodes, edges and prob to sphere around node in row i with param p

% Distance from all nodes to "current" node
thisDistance = pdist2(nodes(i,:), nodes);

thisNodeIdx = thisDistance < p.sphereRadiusOuter & thisDistance > p.sphereRadiusInner;
% rescue inner points that are not connected to center node within inner sphere
innerNodes = thisDistance > p.sphereRadiusInner;
innerEdges = any(ismember(edges, find(innerNodes)),2);
innerConnected = Graph.findConnectedComponents(edges(innerEdges,:));
idx = find(cellfun(@(x)ismember(i,x),innerConnected));
for idx2 = setdiff(1: length(innerConnected), idx)
    thisNodeIdx(innerConnected{idx2}) = true;
end
%fix for not fully connected skeletons
hereConnected = Graph.findConnectedComponents(edges);
idx = cellfun(@(x)ismember(i,x),hereConnected);
if any(idx)
    thisNodeIdx = thisNodeIdx & ismember(1:length(thisDistance),hereConnected{idx});
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