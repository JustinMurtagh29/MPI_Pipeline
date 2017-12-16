function out = patchIntoAgglos(param, agglos, flights)
    % out = patchIntoAgglos(param, agglos, flights)
    %   Patches a set of flight paths (`flights`) into super-agglomerates
    %   (`agglos`).
    %
    %   The `flights` structure must contain a `overlaps` field with a Nx2
    %   matrix indicating the connected agglos. If the entry in the second
    %   column is zero, the corresponding flight is thought to be dangling.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    
    aggloLUT = Superagglos.getSegIds(agglos);
    aggloLUT = Agglo.buildLUT(maxSegId, aggloLUT);
    aggloLUT = cat(1, 0, reshape(aggloLUT, [], 1));
    
    %% grouping agglomerates
    adjEdges = flights.overlaps;
    adjEdges = sort(adjEdges, 2);
    
    % ignore dangling flights
    adjEdges(~adjEdges(:, 1), :) = [];

    adjMat = sparse( ...
        adjEdges(:, 2), adjEdges(:, 1), ...
        true, numel(agglos), numel(agglos));
   [~, aggloComps] = graphconncomp(adjMat, 'Directed', false);

    aggloComps = reshape(aggloComps, [], 1);
    flights.aggloComp = aggloComps(flights.overlaps(:, 1));
    
    %% handle unchanged agglomerates
    changedComps = unique(aggloComps(setdiff(flights.overlaps, 0)));
    unchangedAggloIds = setdiff(1:numel(agglos), flights.overlaps);
    unchangedComps = aggloComps(unchangedAggloIds);
    
    out = struct;
    out.agglos = agglos([]);
    out.agglos(unchangedComps) = agglos(unchangedAggloIds);

    %% patch in flight paths
    for curComp = reshape(changedComps, 1, [])
        curAggloIds = find(aggloComps == curComp);
        curAgglos = agglos(curAggloIds);
        
        curAgglo = struct;
        curAgglo.nodes = cat(1, curAgglos.nodes);
        
        curSegLUT = zeros(maxSegId, 1);
        curNodeIds = find(not(isnan(curAgglo.nodes(:, 4))));
        curSegLUT(curAgglo.nodes(curNodeIds, 4)) = curNodeIds;

        % determine node offset
        curNodeOff = arrayfun(@(a) size(a.nodes, 1), curAgglos);
        curNodeOff = reshape(curNodeOff(1:(end - 1)), [], 1);
        curNodeOff = cumsum(cat(1, 0, curNodeOff));

        curAgglo.edges = cell2mat(arrayfun( ...
            @(ag, off) ag.edges + off, ...
            curAgglos, curNodeOff, 'UniformOutput', false));
        curAgglo.endings = cell2mat(arrayfun( ...
            @(ag, off) ag.endings + off, ...
            curAgglos, curNodeOff, 'UniformOutput', false));

        % patch in flight paths
        curFlightIds = find(flights.aggloComp == curComp);
        curFlightIds = reshape(curFlightIds, 1, []);

        curFlightNodes = cell(numel(curFlightIds), 1);
        curFlightEdges = cell(numel(curFlightIds), 1);
        curAddNodeOff = size(curAgglo.nodes, 1);

        for curIdx = 1:numel(curFlightIds)
            curId = curFlightIds(curIdx);

            curSegIds = cat(2, ...
                flights.segIds{curId}, ...
                flights.neighbours{curId});
            curSegIds = transpose(curSegIds);

            curEndId = aggloLUT(1 + curSegIds);
            curOverlaps = flights.overlaps(curId, :);
            
            if curOverlaps(2)
                % find flight path stretch to extract
                curEndIdx = any(curEndId == curOverlaps(2), 1);
                curEndIdx = 1 + find(curEndIdx(2:end), 1, 'first');
            else
                % flight path is dangling
                curEndIdx = 1 + size(curEndId, 2);
            end
            
            if ~isfield(flights, 'seedNodeIds') ...
                    || ~flights.seedNodeIds(curId)
                % determine start node id
                curStartIdx = any(curEndId == curOverlaps(1), 1);
                curStartIdx = max(find(curStartIdx(1:(curEndIdx - 1)))); %#ok

                curStartNodeId = find( ...
                    curEndId(:, curStartIdx) == curOverlaps(1), 1);
                curStartNodeId = curSegIds(curStartNodeId, curStartIdx);
                curStartNodeId = curSegLUT(curStartNodeId);
            else
                curStartIdx = 0;
                
                % user explicitly passed in ID of seed node
                curStartNodeId = curNodeOff(curAggloIds == curOverlaps(1));
                curStartNodeId = curStartNodeId + flights.seedNodeIds(curId);
            end
            
            if curOverlaps(2)
                curEndNodeId = find( ...
                    curEndId(:, curEndIdx) == curOverlaps(2), 1);
                curEndNodeId = curSegIds(curEndNodeId, curEndIdx);
                curEndNodeId = curSegLUT(curEndNodeId);
            else
                % no end attachment for dangling flights
                curEndNodeId = 0;
            end

            % path in flight nodes
            curNodesToAdd = flights.nodes{curId};
            curNodesToAdd(:, 4) = nan;

            % truncate flights
            curNodesToAdd = curNodesToAdd( ...
                (curStartIdx + 1):(curEndIdx - 1), :);

            % collect new edges
            curNodeCount = size(curNodesToAdd, 1);
            curEdgesToAdd = zeros((curNodeCount - 1) + 2, 2);
            curEdgesToAdd((1 + 1):end, 1) = 1:curNodeCount;
            curEdgesToAdd(1:(end - 1), 2) = 1:curNodeCount;

            curEdgesToAdd = curEdgesToAdd + curAddNodeOff;
            curAddNodeOff = curAddNodeOff + curNodeCount;

            % connect to super-agglomerates
            curEdgesToAdd(1) = curStartNodeId;
            curEdgesToAdd(end) = curEndNodeId;
            
            if ~curEdgesToAdd(end)
                % remove last edge, if dangling
                curEdgesToAdd(end, :) = [];
            end

            curFlightNodes{curIdx} = curNodesToAdd;
            curFlightEdges{curIdx} = curEdgesToAdd;
        end

        curAgglo.nodes = cat(1, curAgglo.nodes, cell2mat(curFlightNodes));
        curAgglo.edges = cat(1, curAgglo.edges, cell2mat(curFlightEdges));
        curAgglo.edges = sort(curAgglo.edges, 2);

        curFlightNodeCount = sum(cellfun( ...
            @(n) size(n, 1), curFlightNodes));
        curAgglo.solvedChiasma = cat( ...
            1, curAgglos.solvedChiasma, ...
            false(curFlightNodeCount, 1));

        % sanity check
        curAdj = sparse( ...
            curAgglo.edges(:, 2), curAgglo.edges(:, 1), ...
            true, size(curAgglo.nodes, 1),size(curAgglo.nodes, 1));
        assert(graphconncomp(curAdj, 'Directed', false) == 1);

        out.agglos(curComp) = curAgglo;
    end
    
    out.agglos = reshape(out.agglos, [], 1);
    out.childIds = reshape(aggloComps, [], 1);
end