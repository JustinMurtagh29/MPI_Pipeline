function superagglosNew = mergeSuperagglosBasedOnFlightPath( ...
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
    
    % Preallocate output
    superagglosNew = struct( ...
        'nodes', cell(numel(eqClassCCfull), 1), ...
        'edges', cell(numel(eqClassCCfull), 1));
    
    % Generate new superagglo
    for i = 1:numel(eqClassCCfull)
        % Concatenate superagglos of this equivalence class
        superagglosNew(i).nodes = cat( ...
            1, superagglos(eqClassCCfull{i}).nodes);
        segmentIds = superagglosNew(i).nodes(:, 4);
        
        % collect and renumber edges
        nrNodes = cumsum(arrayfun(@(x) ...
            size(x.nodes, 1), superagglos(eqClassCCfull{i})));
        nodeOffset = cat(1, 0, nrNodes);
        
        edgeCell = arrayfun( ...
            @(x, y) x.edges + y, ...
            superagglos(eqClassCCfull{i}), nodeOffset(1:end-1), 'uni', 0);
        superagglosNew(i).edges = cat(1, edgeCell{:});
        
        % Find all queries that link to current superagglo (start or end)
        queryIds = find(attachmentsCC == i);
        
        for curFfIdx = queryIds
            % build nodes by adding non-segment ID
            curFfNodes = ff.nodes{curFfIdx};
            curFfNodes(:, 4) = nan;
            
            curFfEdges = Graph.getMST(bsxfun( ...
                @times, curFfNodes(:, 1:3), voxelSize));
            
            % fix node IDs in edge list
            curNodeCount = size(superagglosNew(i).nodes, 1);
            curFfEdges = curFfEdges + curNodeCount;
            
            % build edges connecting query to segment-based agglomerates
            curFfSegIdsHit = cat( ...
                2, ff.segIds{curFfIdx}, ff.neighbours{curFfIdx});
           [~, curFfNodeIdsHit] = ismember(curFfSegIdsHit, segmentIds);
            
            curFfConnEdges = 1:size(curFfNodeIdsHit, 1);
            curFfConnEdges = curFfConnEdges + curNodeCount;
            
            curFfConnEdges = repmat(curFfConnEdges(:), 1, 27);
            curFfConnEdges = cat(2, curFfConnEdges(:), curFfNodeIdsHit(:));
            
            curFfConnEdges = curFfConnEdges(all(curFfConnEdges, 2), :);
            curFfConnEdges = unique(curFfConnEdges, 'rows');
            
            % add to output
            superagglosNew(i).edges = cat( ...
                1, superagglosNew(i).edges, curFfEdges, curFfConnEdges);
            superagglosNew(i).nodes = cat( ...
                1, superagglosNew(i).nodes, curFfNodes);
        end
        
        % sanity checks
        % * edges are correctly sorted
        assert(all(diff(superagglosNew(i,1).edges, 1, 2) > 0));
        
        % * single connected component
        switch size(superagglosNew(i,1).nodes, 1)
            case 0
                assert(false);
            case 1
                assert(isempty(superagglosNew(i,1).edges));
            otherwise
                assert(numel(Graph.findConnectedComponents( ...
                    superagglosNew(i).edges)) == 1);
        end
    end
    
    % make it a column vector
    superagglosNew = reshape(superagglosNew, [], 1);
end
