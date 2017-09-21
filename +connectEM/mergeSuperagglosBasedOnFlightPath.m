function superagglos_new = mergeSuperagglosBasedOnFlightPath( ...
        superagglos, eqClassCCfull, startAgglo, endAgglo, ff)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % configuration
    voxelSize = [11.24, 11.24, 28];
    
    % find connected component for each flight path
    attachments = cellfun(@(x,y)[x y], startAgglo, endAgglo, 'uni', 0);
    ccLookup = Agglo.buildLUT(max(cell2mat(eqClassCCfull)), eqClassCCfull);
    attachmentsCC = cellfun(@(x)unique(ccLookup(x)), attachments);
    
    % Generate new superagglo
    for i=1:length(eqClassCCfull)
        % Concatenate superagglos of this equivalence class
        superagglos_new(i,1).nodes = cat(1, superagglos(eqClassCCfull{i}).nodes);
        
        % collect and renumber edges
        nrNodes = cumsum(arrayfun(@(x)size(x.nodes,1), superagglos(eqClassCCfull{i})));
        nodeOffset = cat(1, 0, nrNodes);
        edgeCell = arrayfun(@(x,y)x.edges+y, superagglos(eqClassCCfull{i}), nodeOffset(1:end-1), 'uni', 0);
        superagglos_new(i,1).edges = cat(1, edgeCell{:});
        
        % Find all queries that link to current superagglo (start or end)
        queryIdx = attachmentsCC == i;
        
        if any(queryIdx)
            % concatenate nodes from flight paths
            % by definition, flight paths nodes have `nan` as segment ID
            newNodes = cat(1, ff.nodes{queryIdx});
            newNodes(:,4) = NaN(size(newNodes,1),1);
            
            % reconstruct flight path edges using MST
            edgeCell = cellfun(@(n) ...
                Graph.getMST(bsxfun(@times, n(:, 1:3), voxelSize)), ...
                ff.nodes(queryIdx), 'uni', 0);
            
            % collect and renumber edges
            nrNodes = cumsum(cellfun(@(x)size(x,1), ff.nodes(queryIdx)));
            nodeOffset = cat(2, 0, nrNodes);
            newEdges = cellfun(@(x,y)x+y, edgeCell, num2cell(nodeOffset(1:end-1)), 'uni', 0);
            newEdges = cat(1, newEdges{:});
            
            % Add in edges that connect current query to current superagglo
            % This means: connect to all segments in neighbourhood
            segIdsThisSuperagglo = superagglos_new(i).nodes(:,4);
            segIdsHitByNewNodes = cat(2, cat(1, ff.segIds{queryIdx}), cat(1, ff.neighbours{queryIdx}));
            [idxNewNodes, idxOldNodes] = ismember(segIdsHitByNewNodes, segIdsThisSuperagglo);
            newNodesBeingAttached = find(any(idxNewNodes,2));
            flatten = @(x)x(:);
            
            % Build segment-to-query edges
            connectingEdges = cat(2, flatten(idxOldNodes(newNodesBeingAttached,:)), repmat(newNodesBeingAttached + size(superagglos_new(i).nodes,1),27,1));
            connectingEdges = unique(connectingEdges, 'rows');
            connectingEdges(any(connectingEdges == 0,2),:) = [];
            
            % Add to output structure
            superagglos_new(i,1).edges = cat(1, superagglos_new(i).edges, newEdges + size(superagglos_new(i).nodes,1), connectingEdges);
            superagglos_new(i,1).nodes = cat(1, superagglos_new(i).nodes, newNodes);
        end
        
        % sanity checks
        % * edges are correctly sorted
        assert(all(diff(superagglos_new(i,1).edges, 1, 2) > 0));
        
        % * single connected component
        switch size(superagglos_new(i,1).nodes, 1)
            case 0
                assert(false);
            case 1
                assert(isempty(superagglos_new(i,1).edges));
            otherwise
                assert(numel(Graph.findConnectedComponents( ...
                    superagglos_new(i).edges)) == 1);
        end
    end
end
