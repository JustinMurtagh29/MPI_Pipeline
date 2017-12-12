function out = patchIntoAgglos(param, axons, flights)
    maxSegId = Seg.Global.getMaxSegId(param);
    
    axonAgglos = Superagglos.getSegIds(axons);
    axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
    axonLUT = cat(1, 0, reshape(axonLUT, [], 1));
    
    %% grouping agglomerates
    adjEdges = flights.overlaps;
    adjEdges = sort(adjEdges, 2);

    adjMat = sparse( ...
        adjEdges(:, 2), adjEdges(:, 1), ...
        true, numel(axons), numel(axons));
   [axonCompCount, axonComps] = ...
        graphconncomp(adjMat, 'Directed', false);

    axonComps = reshape(axonComps, [], 1);
    flights.axonComp = axonComps(flights.overlaps(:, 1));

    %% patch in flight paths
    out = struct;
    out.axons = axons([]);

    for curComp = 1:axonCompCount
        curAxonIds = (axonComps == curComp);
        curAxons = axons(curAxonIds);

        if numel(curAxons) == 1
            out.axons(curComp) = curAxons;
            continue;
        end

        curAxon = struct;
        curAxon.nodes = cat(1, curAxons.nodes);

        curSegLUT = zeros(maxSegId, 1);
        curNodeIds = find(not(isnan(curAxon.nodes(:, 4))));
        curSegLUT(curAxon.nodes(curNodeIds, 4)) = curNodeIds;

        % determine node offset
        curNodeOff = arrayfun(@(a) size(a.nodes, 1), curAxons);
        curNodeOff = cumsum(cat(1, 0, curNodeOff(1:(end - 1))));

        curAxon.edges = cell2mat(arrayfun( ...
            @(ax, off) ax.edges + off, ...
            curAxons, curNodeOff, 'UniformOutput', false));
        curAxon.endings = cell2mat(arrayfun( ...
            @(ax, off) ax.endings + off, ...
            curAxons, curNodeOff, 'UniformOutput', false));

        % patch in flight paths
        curFlightIds = find(flights.axonComp == curComp);
        curFlightIds = reshape(curFlightIds, 1, []);

        curFlightNodes = cell(numel(curFlightIds), 1);
        curFlightEdges = cell(numel(curFlightIds), 1);
        curAddNodeOff = size(curAxon.nodes, 1);

        for curIdx = 1:numel(curFlightIds)
            curId = curFlightIds(curIdx);

            curSegIds = cat(2, ...
                flights.segIds{curId}, ...
                flights.neighbours{curId});
            curSegIds = transpose(curSegIds);

            curEndId = axonLUT(1 + curSegIds);
            curOverlaps = flights.overlaps(curId, :);

            % find flight path stretch to extract
            curEndMask = any(curEndId == curOverlaps(2), 1);
            curEndIdx = 1 + find(curEndMask(2:end), 1, 'first');
            
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
                curStartNodeId = curNodeOff(curAxonIds == curOverlaps(1));
                curStartNodeId = curStartNodeId + flights.seedNodeIds(curId);
            end

            curEndNodeId = find( ...
                curEndId(:, curEndIdx) == curOverlaps(2), 1);
            curEndNodeId = curSegIds(curEndNodeId, curEndIdx);
            curEndNodeId = curSegLUT(curEndNodeId);

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

            curFlightNodes{curIdx} = curNodesToAdd;
            curFlightEdges{curIdx} = curEdgesToAdd;
        end

        curAxon.nodes = cat(1, curAxon.nodes, cell2mat(curFlightNodes));
        curAxon.edges = cat(1, curAxon.edges, cell2mat(curFlightEdges));
        curAxon.edges = sort(curAxon.edges, 2);

        curFlightNodeCount = sum(cellfun( ...
            @(n) size(n, 1), curFlightNodes));
        curAxon.solvedChiasma = cat( ...
            1, curAxons.solvedChiasma, ...
            false(curFlightNodeCount, 1));

        % sanity check
        curAdj = sparse( ...
            curAxon.edges(:, 2), curAxon.edges(:, 1), ...
            true, size(curAxon.nodes, 1),size(curAxon.nodes, 1));
        assert(graphconncomp(curAdj, 'Directed', false) == 1);

        out.axons(curComp) = curAxon;
    end
    
    %{
    otherAxonIds = setdiff(1:numel(allAxons), axonIds);
    out.axons = cat(1, out.axons(:), allAxons(otherAxonIds));

    out.indBigAxons = false(size(out.axons));
    out.indBigAxons(1:axonCompCount) = true;

    % track into which output agglo the input agglos were merged
    out.childIds = nan(numel(allAxons), 1);
    out.childIds(axonIds) = axonComps;
    out.childIds(otherAxonIds) = axonCompCount + (1:numel(otherAxonIds));
    assert(not(any(isnan(out.childIds))));

    % track which flight paths have been used
    out.usedFlightIds = cat(1, flights.id);
    %}
end