function [clusters, centerIds] = ...
        clusterNodes(nodesNm, isIntersection, cutoffDistNm)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    assert(size(nodesNm, 1) == numel(isIntersection));
    nodeIds = find(isIntersection(:));
    
    if isempty(nodeIds)
        clusters = cell(0, 1);
        centerIds = zeros(0, 1);
    elseif isscalar(nodeIds)
        clusters = {nodeIds};
        centerIds = nodeIds;
    else
        clusterIds = clusterdata( ...
            nodesNm(nodeIds, :), ...
            'distance', 'euclidean', ...
            'criterion', 'distance', ...
            'cutoff', cutoffDistNm);

        % build clusters
        clusters = accumarray( ...
            clusterIds, nodeIds, [], @(ids) {ids});

        % find ID of center-most node
       [~, centerIds] = cellfun(@(ids) min(pdist2( ...
            mean(nodesNm(ids, :), 1), nodesNm(ids, :))), clusters);
        centerIds = arrayfun(@(ix, ids) ids{1}(ix), centerIds, clusters);
    end
end