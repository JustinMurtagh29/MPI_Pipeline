function superagglos_new = mergeSuperagglosBasedOnFlightPath(superagglos, eqClassCCfull, startAgglo, endAgglo, ff)

    attachments = cellfun(@(x,y)[x y], startAgglo, endAgglo, 'uni', 0);
    ccLookup = buildLUT(max(cell2mat(eqClassCCfull)), eqClassCCfull);
    attachmentsCC = cellfun(@(x)unique(ccLookup(x)), attachments);
    % Generate new superagglo
    for i=1:length(eqClassCCfull)
        % Concatenate superagglos of this equivalence class
        superagglos_new(i,1).nodes = cat(1, superagglos(eqClassCCfull{i}).nodes);
        nrNodes = cumsum(arrayfun(@(x)size(x.nodes,1), superagglos(eqClassCCfull{i})));
        nodeOffset = cat(1, 0, nrNodes);
        edgeCell = arrayfun(@(x,y)x.edges+y, superagglos(eqClassCCfull{i}), nodeOffset(1:end-1), 'uni', 0);
        superagglos_new(i,1).edges = cat(1, edgeCell{:});
        % Find all queries that link to current superagglo (start or end)
        queryIdx = attachmentsCC == i;
        if any(queryIdx)
            segIdsThisSuperagglo = superagglos_new(i).nodes(:,4);
            % Add in flight path of current query that connects to current eqClass
            newNodes = cat(1, ff.nodes{queryIdx});
            newNodes(:,4) = NaN(size(newNodes,1),1);
            edgeCell = cellfun(@minimalSpanningTree, ff.nodes(queryIdx), 'uni', 0);
            nrNodes = cumsum(cellfun(@(x)size(x,1), ff.nodes(queryIdx)));
            nodeOffset = cat(2, 0, nrNodes);
            newEdges = cellfun(@(x,y)x+y, edgeCell, num2cell(nodeOffset(1:end-1)), 'uni', 0);
            newEdges = cat(1, newEdges{:});
            % Add in edges that connect current query to current superagglo
            segIdsHitByNewNodes = cat(2, cat(1, ff.segIds{queryIdx}), cat(1, ff.neighbours{queryIdx}));
            [idxNewNodes, idxOldNodes] = ismember(segIdsHitByNewNodes, segIdsThisSuperagglo);
            newNodesBeingAttached = find(any(idxNewNodes,2));
            flatten = @(x)x(:);
            connectingEdges = cat(2, flatten(idxOldNodes(newNodesBeingAttached,:)), repmat(newNodesBeingAttached + size(superagglos_new(i).nodes,1),27,1));
            connectingEdges = unique(connectingEdges, 'rows');
            connectingEdges(any(connectingEdges == 0,2),:) = [];
            % Add to output structure
            superagglos_new(i,1).edges = cat(1, superagglos_new(i).edges, newEdges + size(superagglos_new(i).nodes,1), connectingEdges);
            superagglos_new(i,1).nodes = cat(1, superagglos_new(i).nodes, newNodes);
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
        adj(adj > 1000) = 0;
        tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
        [edges(:,1), edges(:,2)] = find(tree);
    end
end

function lut = buildLUT(maxSegId, agglos)
    lut = zeros(maxSegId, 1);
    lut(cell2mat(agglos)) = repelem( ...
        1:numel(agglos), cellfun(@numel, agglos));
end

