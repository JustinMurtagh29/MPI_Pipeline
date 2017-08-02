function [pos, dir, queryIdx] = detectChiasmataPostSingleNodeLabelSubSub(node,edges,prob,p,cc,centerOfCC,pos,dir,queryIdx,i)
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
        queryTemp=find(ismember(nodes, thisNodes(closestNode, :),'rows'));
        queryIdx{i}(end + 1) = max([-1;queryTemp(:)]);
    end
end
function dummy()
end