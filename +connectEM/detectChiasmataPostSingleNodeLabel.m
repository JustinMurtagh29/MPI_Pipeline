function [output, queryIdx] = detectChiasmataPostSingleNodeLabel(edges, isIntersection, nrExits, nodes, p, nodesV, prob)
temp.edges = edges;
cc = findCCaccordingToGraph(temp, find(isIntersection)); %must be from manuelCode repo
[~, centerOfCC] = cellfun(@(x)min(pdist2(bsxfun(@minus, nodes(x,:), mean(nodes(x,:),1)), [0 0 0])), cc);

% Find out where to query for each CC
queryIdx = cell(length(cc),1);
pos = cell(length(cc),1);
dir = cell(length(cc),1);
for i=1:length(cc)
    i
    
    [thisNodes, thisEdges, thisProb] = connectEM.detectChiasmataPruneToSphere(nodes, edges, prob, p, cc{i}(centerOfCC(i)));
    C = Graph.findConnectedComponents(thisEdges);
    goodcomps = find(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodes(i,:))) > 4000, C));
    assert(size(goodcomps,2) == 1); % make sure that the next line works
    for idx = goodcomps'
        thisNodesCC = inf(size(thisNodes));
        thisNodesCC(C{idx},:) = thisNodes(C{idx}, :);
        [~, closestNode] = min(pdist2(thisNodesCC, nodes(cc{i}(centerOfCC(i)),:)));
        bestEdge = find(any(thisEdges == closestNode, 2), 1);
        descale = @(x)bsxfun(@times, x, 1./p.voxelSize);
        pos{i}(end + 1, :) = descale(thisNodes(closestNode,:));
        dir{i}(end + 1, :) = bsxfun(@minus, pos{i}(end, :), descale(thisNodes(setdiff(thisEdges(bestEdge,:), closestNode),:)));
        queryIdx{i}(end + 1) =find(ismember(nodes, thisNodes(closestNode, :),'rows'));
        
    end
    
end

% Create an output structure
output.nodes = nodesV;
output.edges = edges;
output.prob = [];
output.isIntersection = isIntersection;
output.nrExits = nrExits;
output.ccNodeIdx = cc;
output.ccCenterIdx = cellfun(@(x,y)x(y), cc, num2cell(centerOfCC));
output.queryIdx = queryIdx;
output.position = pos;
output.direction = dir;

