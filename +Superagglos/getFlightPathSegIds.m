function flights = getFlightPathSegIds(param, sagglos, nhood)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('nhood', 'var') || isempty(nhood)
        nhood = 1;
    else
        assert(mod(nhood, 2) == 1);
    end
    
    %% reconstruct flights
    aggloCount = numel(sagglos);
    flights = cell(aggloCount, 1);
    flights(:) = {zeros(0, 6)};

    for curIdx = 1:aggloCount
        curSagglos = sagglos(curIdx);

        % work around empty edges
        curSagglos.edges = reshape(curSagglos.edges, [], 2);

        % extract flight nodes and edges
        curNodeIds = find(isnan(curSagglos.nodes(:, 4)));
       [~, curEdges] = ismember(curSagglos.edges, curNodeIds);
        curEdges = curEdges(all(curEdges, 2), :);
        curEdges = sort(curEdges, 2);

        % group into paths
        curAdj = sparse( ...
            curEdges(:, 2), curEdges(:, 1), ...
            true, numel(curNodeIds), numel(curNodeIds));
       [~, curLUT] = graphconncomp(curAdj, 'Directed', false);

        flights{curIdx} = horzcat( ...
           repelem(curIdx, numel(curNodeIds), 1), curNodeIds(:), ...
           curLUT(:), curSagglos.nodes(curNodeIds, 1:3));
    end

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