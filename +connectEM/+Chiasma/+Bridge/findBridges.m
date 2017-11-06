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
    
    %%
    curIdA = exitNodeIds(1);
    for curIdB = exitNodeIds(2:end)
        curIdsAB = horzcat(curIdA, curIdB);
        curIdsCD = setdiff(exitNodeIds, curIdsAB);
        
        % add loops
        curEdges = axon.edges;
        curEdges = cat(1, curEdges, [curIdA, curIdB; curIdsCD]);
        assert(issorted(curEdges, 2));
        
        curAdj = sparse( ...
            curEdges(:, 2), curEdges(:, 1), ...
            true, nodeCount, nodeCount);
        assert(graphconncomp(curAdj, 'Directed', false) == 1);

        % search for loops
        for curEdgeId = 1:size(curEdges, 1)
            curAdj = curEdges;
            curAdj(curEdgeId, :) = [];
            assert(issorted(curAdj, 2));
        
            curAdj = sparse( ...
                curAdj(:, 2), curAdj(:, 1), ...
                true, nodeCount, nodeCount);
           [~, curComps] = graphconncomp(curAdj, 'Directed', false);
           
            curCompAB = unique(curComps(curIdsAB));
            curCompCD = unique(curComps(curIdsCD));
            
            assert(isscalar(curCompAB));
            assert(isscalar(curCompCD));
            
            if curCompAB ~= curCompCD
                fprintf('Found bridge:\n');
                fprintf('(%d, %d) and (%d, %d)\n', curIdsAB, curIdsCD);
                fprintf('Edge: (%d, %d)\n', curEdges(curEdgeId, :));
            end
        end
    end
end