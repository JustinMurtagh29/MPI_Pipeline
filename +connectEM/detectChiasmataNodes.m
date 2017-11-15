function [nrExits, pos, dir, queryIdx] = ...
        detectChiasmataNodes(p, nodes, edges, nodeIdx)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
   [~, thisEdges, thisNodeIds, thisDist] = ...
        connectEM.detectChiasmataPruneToSphere(p, nodes, edges, nodeIdx);
   [clusters, clusterIds] = ...
        Graph.findConnectedComponents(thisEdges, false);
   
    % restrict to true exits
    clusters = clusters(accumarray( ...
        clusterIds(:), thisDist, [], @max) > p.minNodeDist);
    nrExits = numel(clusters);
    clear lut;
    
    %% do query generation, if desired
    if nargout < 2; return; end
    
    pos = nan(nrExits, 3);
    dir = nan(nrExits, 3);
    queryIdx = nan(1, nrExits);
    
    % NOTE(amotta): At least one end of the query edges must reside within
    % the remains of `detectChiasmataPruneToSphere`.
    mask = any(ismember(edges, thisNodeIds), 2);
    edges = edges(mask, :);
    clear mask;
    
    % NOTE(amotta): Both ends of the query edge must be less than 9 µm from
    % the current node-of-interest.
    mask = unique(edges);
    mask = mask(pdist(nodes(mask, :), nodes(nodeIdx, :)) < 9000);
    mask = all(ismember(edges, mask), 2);
    edges = edges(mask, :);
    clear mask;
    
    descale = @(x) bsxfun(@times, x, 1 ./ p.raw.voxelSize);
    
    for curIdx = 1:numel(clusters)
        curNodeIds = thisNodeIds(clusters{curIdx});
        
        curMask = ismember(edges, curNodeIds);
        curMask = find(xor(curMask(:, 1), curMask(:, 2)), 1);
        curNodeId = intersect(edges(curMask, :), curNodeIds);
        
        pos(curIdx, :) = descale(nodes(curNodeId,:));
        dir(curIdx, :) = pos(curIdx, :) - descale(nodes(nodeIdx, :));
        queryIdx(curIdx) = curNodeId;
    end
    
    % sanity checks
    assert(size(pos, 1) == nrExits);
    assert(size(dir, 1) == nrExits);
    assert(size(queryIdx, 2) == nrExits);
end