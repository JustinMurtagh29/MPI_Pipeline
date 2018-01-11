function [flights, flightsMeta] = ...
        getFlightPathSegIds(param, sagglos, nhood)
    % [flights, flightsMeta] = getFlightPathSegIds(param, sagglos, nhood)
    %   This function takes a bunch of super-agglomerates `sagglos`, splits
    %   them in connected segment- and flight-based parts, and looks up the
    %   segment IDs for all flight nodes. Optionally, the segment IDs are
    %   being looked up for a whole neighborhood around the node.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('nhood', 'var') || isempty(nhood)
        nhood = 1;
    else
        assert(mod(nhood, 2) == 1);
    end
    
    %% reconstruct flights
    aggloCount = numel(sagglos);
    flights = repelem({zeros(0, 6)}, aggloCount, 1);
    flightsMeta = repelem({zeros(0, 3)}, aggloCount, 1);

    for curIdx = 1:aggloCount
        curSagglo = sagglos(curIdx);
        curNodeCount = size(curSagglo.nodes, 1);

        % work around empty edges
        curSagglo.edges = reshape(curSagglo.edges, [], 2);
        
        % separate flight paths from agglomerates
        curFlightNodes = find(isnan(curSagglo.nodes(:, 4)));
        curInterEdges = ismember(curSagglo.edges, curFlightNodes);
        curInterEdges = xor(curInterEdges(:, 1), curInterEdges(:, 2));
        
        curIntraEdges = curSagglo.edges(~curInterEdges, :);
        curInterEdges = curSagglo.edges( curInterEdges, :);

        % build agglomerates
        curAdj = sparse( ...
            curIntraEdges(:, 2), curIntraEdges(:, 1), ...
            true, curNodeCount, curNodeCount);
       [~, curLUT] = graphconncomp(curAdj, 'Directed', false);
        assert(all(curLUT > 0));
        
        % extract
        % * nodes for each flight path
        % * number of attachments per flight
        curInterEdges = curLUT(curInterEdges);
        curInterEdges = sort(curInterEdges, 2);
        curInterEdges = unique(curInterEdges, 'rows');
        
       [curFlightComps, ~, curLUT] = unique(curLUT(curFlightNodes));
       [~, curFlightAttachments] = ismember(curInterEdges, curFlightComps);
        curFlightAttachments = accumarray( ...
            max(curFlightAttachments, [], 2), ...
            1, [numel(curFlightComps), 1], [], 0);
        
        flights{curIdx} = horzcat( ...
            repelem(curIdx, numel(curFlightNodes), 1), curFlightNodes(:), ...
            curLUT(:), curSagglo.nodes(curFlightNodes, 1:3));
        flightsMeta{curIdx} = horzcat( ...
            repelem(curIdx, numel(curFlightComps), 1), ...
            reshape(1:numel(curFlightComps), [], 1), ...
            reshape(curFlightAttachments, [], 1));
    end
    
    flightsMeta = array2table( ...
        cell2mat(flightsMeta), 'VariableNames', ...
       {'aggloId', 'flightId', 'numAttachments'});

    flights = cell2mat(flights);
    flights = transpose(flights);

    %% add 26-neighbours
    if nhood > 1
        nhoodVec = (nhood - 1) / 2;
        nhoodVec = (-nhoodVec):nhoodVec;
        
        posOff = cell(3, 1);
       [posOff{:}] = ndgrid(nhoodVec, nhoodVec, nhoodVec);
       
        posOff = cat(4, posOff{:});
        posOff = reshape(posOff, [], 3);
        posOff = transpose(posOff);

        posOff = cat(1, zeros( ...
            size(flights, 1) - size(posOff, 1), ...
            size(posOff, 2)), posOff);
        
        flights = reshape(flights, 6, 1, []);
        flights = bsxfun(@plus, flights, posOff);
        flights = reshape(flights, 6, []);
        flights = transpose(flights);
    end
    
    %% convert into table
    nodeCoord = flights(:, (end - 2):end);
    flights = flights(:, 1:(end - 3));
    
    flights = array2table( ...
        flights, 'VariableNames', ...
       {'aggloId', 'nodeId', 'flightId'});

    %% look up segment IDs
    flights.segId = Seg.Global.getSegIds( ...
        param, nodeCoord, [1024, 1024, 1024]);
end