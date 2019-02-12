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
        
        curSegNodeMask = not(isnan(curSagglo.nodes(:, 4)));
        
        % Segment-based parts
        curNodeIds = find(curSegNodeMask);
       [~, curEdges] = ismember(curSagglo.edges, curNodeIds);
        curEdges = curEdges(all(curEdges, 2), :);
        
        agglos{curIdx} = doIt(curSagglo, curNodeIds, curEdges);
        
        % Flight-based parts
        curNodeIds = not(all(curSegNodeMask(curSagglo.edges), 2));
        curNodeIds = unique(curSagglo.edges(curNodeIds, :));
       [~, curEdges] = ismember(curSagglo.edges, curNodeIds);
        curEdges = curEdges(all(curEdges, 2), :);
        
        flights{curIdx} = doIt(curSagglo, curNodeIds, curEdges);
    end
end

function skels = doIt(sagglo, nodeIds, relEdges)
    curGraph = graph(relEdges(:, 1), relEdges(:, 2), [], numel(nodeIds));
    curComps = curGraph.conncomp('OutputForm', 'cell');
    
    skels = struct('nodes', {}, 'edges', {});
    skels = repelem(skels, numel(nodeIds), 1);
    
    for curId = 1:numel(curComps)
        curRelIds = curComps{curId};
        
        curNodes = nodeIds(curRelIds);
        curNodes = sagglo.nodes(curNodes, :);
        skels(curId).nodes = curNodes;
        
       [~, curEdges] = ismember(relEdges, curRelIds);
        curEdges = curEdges(all(curEdges, 2), :);
        skels(curId).edges = curEdges;
    end
end
