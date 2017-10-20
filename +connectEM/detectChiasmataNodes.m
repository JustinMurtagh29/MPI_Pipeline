function [nrExits, pos, dir, queryIdx] = ...
        detectChiasmataNodes(nodes, edges, prob, p, nodeIdx)
    [thisNodes, thisEdges] = connectEM.detectChiasmataPruneToSphere( ...
            nodes, edges, prob, p, nodeIdx);
    C = Graph.findConnectedComponents(thisEdges, false);
    C = C(cellfun(@(idx) max(pdist2( ...
        thisNodes(idx, :), nodes(nodeIdx, :))) > p.minDistNode, C));
    nrExits = numel(C);
    
    %% do query generation, if desired
    if nargout < 2; return; end
    
    pos = nan(nrExits, 3);
    dir = nan(nrExits, 3);
    queryIdx = nan(1, nrExits);
    
    nodesMask = find(ismember(nodes, thisNodes, 'rows'));
    
    for idx = 1:numel(C)
        nodesMask2 = find(ismember(nodes, thisNodes(C{idx}, :), 'rows'));
        transitioningEdge = find(any(ismember(edges,nodesMask2),2)& ... %one node of the edge should be a member of nodesMask2
            any(~ismember(edges,nodesMask),2)& ... %the other shouldn't be a member of nodesMask or nodesMask2 (which is a subset of nodesMask)
            all(ismember(edges,find(pdist2(nodes,nodes(nodeIdx,:))<(p.sphereRadiusOuter))),2),1); % and we don't want an edge from the outer cutting
        descale = @(x)bsxfun(@times, x, 1./p.voxelSize);
        queryTemp=intersect(edges(transitioningEdge,:), nodesMask2);
        pos(idx, :) = descale(nodes(queryTemp,:));
        dir(idx, :) = bsxfun(@minus, pos(idx, :), descale(nodes(nodeIdx,:)));
        queryIdx(idx) = max([-1;queryTemp(:)]);
    end
    
    % sanity checks
    assert(size(pos, 1) == nrExits);
    assert(size(dir, 1) == nrExits);
    assert(size(queryIdx, 2) == nrExits);
end