function superagglos_new = mergeSuperagglosBasedOnFlightPath(superagglos, eqClassCCfull, startAgglo, endAgglo, ff)

    attachments = cellfun(@(x,y)[x y], startAgglo, endAgglo, 'uni', 0);
    % Generate new superagglo
    for i=1:length(eqClassCCfull)
        % Concatenate superagglos of this equivalence class
        superagglos_new(i).nodes = cat(1, superagglos(eqClassCCfull{i}).nodes);
        nrNodes = arrayfun(@(x)size(x.nodes,1), superagglos(eqClassCCfull{i}));
        nodeOffset = cat(2, 0, nrNodes(1:end-1));
        edgeCell = arrayfun(@(x,y)x.edges+y, superagglos(eqClassCCfull{i}), nodeOffset, 'uni', 0);
        superagglos_new(i).edges = cat(1, edgeCell{:});
        % Find all queries that link to current superagglo (start or end) and all segIds of superagglo
        queryIdx = find(cellfun(@(x)any(ismember(x, eqClassCCfull{i})), attachments));
        
        segIdsThisSuperagglo = superagglos_new(i).nodes(:,4);
        for j=1:length(queryIdx)
            % Add in flight path of current query that connects to current eqClass
            newNodes = ff.nodes{queryIdx(j)};
            newNodes(:,4) = NaN(size(newNodes,1),1);
            newEdges = minimalSpanningTree(newNodes(:,1:3));
            segIdsHitByNewNodes = cat(2, ff.segIds{queryIdx(j)}, ff.neighbours{queryIdx(j)});
            % Add in edges that connect current query to current superagglo
            [idxNewNodes, idxOldNodes] = ismember(segIdsHitByNewNodes, segIdsThisSuperagglo);
            newNodesBeingAttached = find(any(idxNewNodes,2));
            flatten = @(x)x(:);
            connectingEdges = cat(2, flatten(idxOldNodes(newNodesBeingAttached,:)), repmat(newNodesBeingAttached + size(superagglos_new(i).nodes,1),27,1));
            connectingEdges(any(connectingEdges == 0,2),:) = [];
            % Add to output structure
            superagglos_new(i).edges = cat(1, superagglos_new.edges, cat(1, newEdges + size(superagglos_new(i).nodes,1), connectingEdges));
            superagglos_new(i).nodes = cat(1, superagglos_new.nodes, newNodes);
        end
    end
end

function edges = minimalSpanningTree(com)
    com = bsxfun(@times, com, [11.24 11.24 28]);
    if size(com,1) < 2
        edges = [];
    else
        % Minimal spanning tree, probably inefficent
        adj = squareform(pdist(com));
        tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
        [edges(:,1), edges(:,2)] = find(tree);
    end
end

