function [agglos, flights] = splitIntoAgglosAndFlights(sagglos)
    % [agglos, flights] = splitIntoAgglosAndFlights(sagglos)
    %   Splits a super-agglomerate into the underlying (segment-based)
    %   agglomerates and flight path parts.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % prepare output
    template = struct('nodes', {}, 'edges', {});
    template = {reshape(template, [], 1)};
    
    agglos = repmat(template, size(sagglos));
    flights = repmat(template, size(sagglos));
    
    for curIdx = 1:numel(sagglos)
        curSagglo = sagglos(curIdx);
        assert(issorted(curSagglo.edges, 2));

        % separate flight paths from agglomerates
        curIntraEdges = find(isnan(curSagglo.nodes(:, 4)));
        curIntraEdges = ismember(curSagglo.edges, curIntraEdges);
        curIntraEdges = curSagglo.edges(sum(curIntraEdges, 2) ~= 1, :);

       [curNodeIds, ~, curIntraEdges] = unique(curIntraEdges);
        curIntraEdges = reshape(curIntraEdges, [], 2);

        curNodeCount = numel(curNodeIds);
        curNodes = curSagglo.nodes(curNodeIds, :);

        % find components
        curAdj = sparse( ...
            curIntraEdges(:, 2), curIntraEdges(:, 1), ...
            true, curNodeCount, curNodeCount);
       [curCompCount, curLUT] = ...
            graphconncomp(curAdj, 'Directed', false);
        assert(all(curLUT > 0));

        % extract components
        for curCompIdx = 1:curCompCount
            curCompNodeIds = find(curLUT == curCompIdx);
            curCompNodes = curNodes(curCompNodeIds, :);

            % edges
           [~, curCompEdges] = ismember(curIntraEdges, curCompNodeIds);
            curCompEdges(~all(curCompEdges, 2), :) = [];
            
            curCompSagglo = struct( ...
                'nodes', curCompNodes, ...
                'edges', curCompEdges);
            
            if any(isnan(curNodes(curCompNodeIds, 4)))
                flights{curIdx}(end + 1, 1) = curCompSagglo;
            else
                agglos{curIdx}(end + 1, 1) = curCompSagglo;
            end
        end
    end
end