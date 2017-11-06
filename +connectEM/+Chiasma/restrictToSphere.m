function [postAgglo, preNodeIds, preNodeDistNmSq] = ...
        restrictToSphere(param, agglo, nodeIdx, radiusNm)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~isfield(agglo, 'nodesNm')
        voxelSize = param.raw.voxelSize;
        agglo.nodesNm = agglo.nodes(:, 1:3) .* voxelSize;
    end
    
    % restrict spatially
    nodePosNm = agglo.nodesNm(nodeIdx, :);
    preNodeDistNmSq = pdist2(agglo.nodesNm, nodePosNm, 'squaredeuclidean');
    preNodeIds = find(preNodeDistNmSq(:) < radiusNm * radiusNm);
    
    % restrict to center component
   [~, postEdges] = ismember(agglo.edges, preNodeIds);
    postEdges = postEdges(all(postEdges, 2), :);
    
    curComps = sparse( ...
        postEdges(:, 2), postEdges(:, 1), ...
        true, numel(preNodeIds), numel(preNodeIds));
   [~, curComps] = graphconncomp(curComps, 'Directed', false);
    
    curCompId = curComps(preNodeIds == nodeIdx);
    preNodeIds = preNodeIds(curComps == curCompId);
    
    % build output
    postAgglo = struct;
    postAgglo.nodes = agglo.nodes(preNodeIds, :);
    postAgglo.nodesNm = agglo.nodesNm(preNodeIds, :);
    
   [~, postAgglo.edges] = ismember(agglo.edges, preNodeIds);
    postAgglo.edges = postAgglo.edges(all(postAgglo.edges, 2), :);
end