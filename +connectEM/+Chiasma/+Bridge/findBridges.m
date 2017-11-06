function findBridges(param, chiasmaParam, axon, chiasmata, chiasmaId)
    curNodeId = chiasmata.ccCenterIdx(chiasmaId);
    curExitNodeIds = chiasmata.queryIdx{chiasmaId};
    
   [~, curInnerNodeIds] = connectEM.Chiasma.restrictToSphere( ...
       param, axon, curNodeId, chiasmaParam.sphereRadiusInner);
    
    % add immediate neighbours as well
    curNodeIds = any(ismember(axon.edges, curInnerNodeIds), 2);
    curNodeIds = unique(axon.edges(curNodeIds, :));
    curNodeCount = numel(curNodeIds);
    
    % restrict
    axon.nodes = axon.nodes(curNodeIds, :);
   [~, axon.edges] = ismember(axon.edges, curNodeIds);
    axon.edges(~all(axon.edges, 2), :) = [];
    
   [~, curExitNodeIds] = ismember(curExitNodeIds, curNodeIds);
    curExitNodeIds = sort(reshape(curExitNodeIds, 1, 4));
    assert(all(curExitNodeIds));
    
    %%
    curIdA = curExitNodeIds(1);
    for curIdB = curExitNodeIds(2:end)
        curIdsAB = horzcat(curIdA, curIdB);
        curIdsCD = setdiff(curExitNodeIds, curIdsAB);
        
        % add loops
        curEdges = axon.edges;
        curEdges = cat(1, curEdges, [curIdA, curIdB; curIdsCD]);
        assert(issorted(curEdges, 2));
        
        curAdj = sparse( ...
            curEdges(:, 2), curEdges(:, 1), ...
            true, curNodeCount, curNodeCount);
        assert(graphconncomp(curAdj, 'Directed', false) == 1);

        % search for loops
        for curEdgeId = 1:size(curEdges, 1)
            curAdj = curEdges;
            curAdj(curEdgeId, :) = [];
            assert(issorted(curAdj, 2));
        
            curAdj = sparse( ...
                curAdj(:, 2), curAdj(:, 1), ...
                true, curNodeCount, curNodeCount);
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