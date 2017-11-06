function findBridges(param, chiasmaParam, axon, chiasmata, chiasmaId)
    nodeId = chiasmata.ccCenterIdx(chiasmaId);
    exitNodeIds = chiasmata.queryIdx{chiasmaId};
    
   [~, innerNodeIds] = connectEM.Chiasma.restrictToSphere( ...
       param, axon, nodeId, chiasmaParam.sphereRadiusInner);
    
    % add immediate neighbours as well
    nodeIds = any(ismember(axon.edges, innerNodeIds), 2);
    nodeIds = unique(axon.edges(nodeIds, :));
    nodeCount = numel(nodeIds);
    
    % restrict
    axon.nodes = axon.nodes(nodeIds, :);
   [~, axon.edges] = ismember(axon.edges, nodeIds);
    axon.edges(~all(axon.edges, 2), :) = [];
    
   [~, exitNodeIds] = ismember(exitNodeIds, nodeIds);
    exitNodeIds = sort(reshape(exitNodeIds, 1, 4));
    assert(all(exitNodeIds));
    
    %% find bridges
    edges = axon.edges;
    assert(issorted(edges, 2));
    
    adjMat = sparse( ...
        edges(:, 2), edges(:, 1), ...
        true, nodeCount, nodeCount);
    assert(graphconncomp(adjMat, 'Directed', false) == 1);
    clear adjMat;
    
    % search for bridges
    for curEdgeId = 1:size(edges, 1)
        curEdges = edges;
        curEdges(curEdgeId, :) = [];
        
        curAdjMat = sparse( ...
            curEdges(:, 2), curEdges(:, 1), ...
            true, nodeCount, nodeCount);
        [~, curComps] = graphconncomp(curAdjMat, 'Directed', false);
        
        % check for 2-2 partition
        for curIdB = exitNodeIds(2:end)
            curIdsAB = horzcat(exitNodeIds(1), curIdB);
            curIdsCD = setdiff(exitNodeIds, curIdsAB);

            curCompAB = unique(curComps(curIdsAB));
            if ~isscalar(curCompAB); continue; end
            
            curCompCD = unique(curComps(curIdsCD));
            if ~isscalar(curCompCD); continue; end
            
            if curCompAB ~= curCompCD
                fprintf('Found bridge:\n');
                fprintf('(%d, %d) and (%d, %d)\n', curIdsAB, curIdsCD);
                fprintf('Edge: (%d, %d)\n', edges(curEdgeId, :));
            end
        end
    end
end