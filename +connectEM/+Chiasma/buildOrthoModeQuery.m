function skel = buildOrthoModeQueryFromDetection(param, agglo, nodeIdx)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % main logic
    agglo.nodesNm = agglo.nodes(:, 1:3) .* param.raw.voxelSize;
    centerNodeNm = agglo.nodesNm(nodeIdx, :);
    
    % generate fake edge probabilities
    % These are not actually of any importance.
    edgeProb = ones(size(agglo.edges, 1), 1);
    
    % cut out sphere
   [cutNodes, cutEdges, ~, cutNodeIds] = ...
       connectEM.detectChiasmataPruneToSphere( ...
            agglo.nodesNm, agglo.edges, edgeProb, param, nodeIdx);    
	cutComps = Graph.findConnectedComponents(cutEdges, false);
    
    % drop tiny components
    cutComps = cutComps(cellfun(@(idx) max(pdist2( ...
        cutNodes(idx, :), centerNodeNm)) > param.minNodeDist, cutComps));
    cutCompCount = numel(cutComps);
    
    % generate color
    cutColors = distinguishable_colors(cutCompCount, [1 1 1; 0 0 0]);
    cutColors = cat(2, cutColors, ones(cutCompCount, 1));
    
    % build skeleton
    skel = skeleton();
    
    for curIdx = 1:cutCompCount
        curNodeIds = cutComps{curIdx};
        
        % relabel edges
       [~, curEdges] = ismember(cutEdges, curNodeIds);
        curEdges = curEdges(all(curEdges, 2), :);
        
        % get nodes
        curNodeIds = cutNodeIds(curNodeIds);
        curNodes = agglo.nodes(curNodeIds, 1:3);
        
        % get color
        curColor = cutColors(curIdx, :);
        
       [~, curCommentIdx] = min(pdist2( ...
           agglo.nodesNm(curNodeIds, :), centerNodeNm));
        curComments = repelem({''}, numel(curNodeIds), 1);
        curComments{curCommentIdx} = sprintf('Exit %d', curIdx);
        
        skel = skel.addTree( ...
            sprintf('Branch %d', curIdx), ...
            curNodes, curEdges, curColor, [], curComments);
    end
    
    skel = Skeleton.setParams4Pipeline(skel, param);
end