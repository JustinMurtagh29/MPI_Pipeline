function [isBridge, edges] = findBridges( ...
        param, chiasmaParam, axon, chiasmata, chiasmaId)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    nodeId = chiasmata.ccCenterIdx(chiasmaId);
    exitNodeIds = chiasmata.queryIdx{chiasmaId};
    
    % restrict to outer sphere
   [axon, outerNodeIds] = connectEM.Chiasma.restrictToSphere( ...
       param, axon, nodeId, chiasmaParam.sphereRadiusOuter);
   
    nodeId = find(outerNodeIds == nodeId);
   [~, exitNodeIds] = ismember(exitNodeIds, outerNodeIds);
   
    % find nodes which are cut out
   [~, innerNodeIds] = connectEM.Chiasma.restrictToSphere( ...
       param, axon, nodeId, chiasmaParam.sphereRadiusInner); %#ok
    
    % add immediate neighbours as well
    nodeIds = any(ismember(axon.edges, innerNodeIds), 2);
    nodeIds = unique(axon.edges(nodeIds, :));
    nodeCount = numel(nodeIds);
    
    % restrict
    axon.nodes = axon.nodes(nodeIds, :);
   [~, axon.edges] = ismember(axon.edges, nodeIds);
    axon.edges(~all(axon.edges, 2), :) = [];
    
   [~, exitNodeIds] = ismember(exitNodeIds, outerNodeIds(nodeIds));
    exitNodeIds = reshape(exitNodeIds, 1, 4);
    assert(all(exitNodeIds));
    
    %% find bridges
    edges = axon.edges;
    assert(issorted(edges, 2));
    
    edgeCount = size(edges, 1);
    isBridge = false(edgeCount, 3);
    
    adjMat = sparse( ...
        edges(:, 2), edges(:, 1), ...
        true, nodeCount, nodeCount);
    assert(graphconncomp(adjMat, 'Directed', false) == 1);
    clear adjMat;
    
    % search for bridges
    for curEdgeId = 1:edgeCount
        curEdges = edges;
        curEdges(curEdgeId, :) = [];
        
        curAdjMat = sparse( ...
            curEdges(:, 2), curEdges(:, 1), ...
            true, nodeCount, nodeCount);
        [~, curComps] = graphconncomp(curAdjMat, 'Directed', false);
        
        for curIdxB = 2:4
            curIdsAB = sort(exitNodeIds([1, curIdxB]));
            curIdsCD = setdiff(exitNodeIds, curIdsAB);

            curCompAB = unique(curComps(curIdsAB));
            if ~isscalar(curCompAB); continue; end
            
            curCompCD = unique(curComps(curIdsCD));
            if ~isscalar(curCompCD); continue; end
            
            % check for 2-2 partition
            if curCompAB == curCompCD; continue; end
            isBridge(curEdgeId, curIdxB - 1) = true;
        end
    end
    
    % globalize edges
    edges = outerNodeIds(nodeIds(edges));
end