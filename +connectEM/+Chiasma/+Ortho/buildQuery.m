function [skel, query] = buildQuery(param, agglo, nodeIdx)
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
    activeNodeId = nan;
    cutExitNodeIds = nan(cutCompCount, 1);
    
    for curIdx = 1:cutCompCount
        curNodeIds = cutComps{curIdx};
        
        % relabel edges
       [~, curEdges] = ismember(cutEdges, curNodeIds);
        curEdges = curEdges(all(curEdges, 2), :);
        
        % get nodes
        curNodeIds = cutNodeIds(curNodeIds);
        curNodes = agglo.nodes(curNodeIds, 1:3);
        
        % find exit node
       [~, curExitNodeId] = min(pdist2( ...
           agglo.nodesNm(curNodeIds, :), centerNodeNm));
        
        % retain exit nodes and mark active node
        cutExitNodeIds(curIdx) = curNodeIds(curExitNodeId);
        if curIdx == 1; activeNodeId = curExitNodeId; end
        
        % get color
        curColor = cutColors(curIdx, :);
       
        curComments = repelem({''}, numel(curNodeIds), 1);
        curComments{curExitNodeId} = sprintf('Exit %d', curIdx);
        
        skel = skel.addTree( ...
            sprintf('Branch %d', curIdx), ...
            curNodes, curEdges, curColor, [], curComments);
    end
    
    % set dataset
    skel = Skeleton.setParams4Pipeline(skel, param);
    
    % set viewing positions
    skel.parameters.editPosition = struct( ...
        'x', agglo.nodes(cutExitNodeIds(1), 1), ...
        'y', agglo.nodes(cutExitNodeIds(1), 2), ...
        'z', agglo.nodes(cutExitNodeIds(1), 3));
    skel.parameters.zoomLevel.zoom = 0.53;
    
    % set active node
    skel.parameters.activeNode.id = activeNodeId;
    
    % build `query`
    query = struct;
    query.centerNodeId = nodeIdx;
    query.exitNodeIds = cutExitNodeIds;
end