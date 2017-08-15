function [pos, dir, queryIdx] = detectChiasmataPostSingleNodeLabelSubSub(nodes,edges,prob,p,cc,centerOfCC,pos,dir,queryIdx,i)
    [thisNodes, thisEdges, thisProb] = connectEM.detectChiasmataPruneToSphere(nodes, edges, prob, p, cc{i}(centerOfCC(i)));
    nodesMask = find(ismember(nodes, thisNodes, 'rows'));
    C = Graph.findConnectedComponents(thisEdges, false);
    goodcomps = find(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodes(cc{i}(centerOfCC(i)),:))) > 3000, C));
    assert(size(goodcomps,2) == 1); % make sure that the next line works

        
        

    for idx = goodcomps'
        nodesMask2 = find(ismember(nodes, thisNodes(C{idx}, :), 'rows'));
        transitioningEdge = find(any(ismember(edges,nodesMask2),2)& ... %one node of the edge should be a member of nodesMask2
            any(~ismember(edges,nodesMask),2)& ... %the other shouldn't be a member of nodesMask or nodesMask2 (which is a subset of nodesMask)
            all(ismember(edges,find(pdist2(nodes,nodes(cc{i}(centerOfCC(i)),:))<9000)),2),1); % and we don't want an edge from the outer cutting
        descale = @(x)bsxfun(@times, x, 1./p.voxelSize);
        queryTemp=intersect(edges(transitioningEdge,:), nodesMask2);
        pos{i}(end + 1, :) = descale(nodes(queryTemp,:));
        dir{i}(end + 1, :) = bsxfun(@minus, pos{i}(end, :), descale(nodes(cc{i}(centerOfCC(i)),:)));
        queryIdx{i}(end + 1) = max([-1;queryTemp(:)]);
    end
end
function dummy()
end