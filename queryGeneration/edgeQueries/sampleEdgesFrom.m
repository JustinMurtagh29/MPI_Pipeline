function edges = sampleEdgesFrom(graph, segIds, skipEdges)
    % edges = sampleEdgesFrom(graph, segIds, skipEdges)
    %   Randomly samples edges from the whole data set, which
    %   can then be used to build edge queries.
    %
    % graph
    %   Graph structure with the 'edges' property
    %
    % segIds
    %   List of global segment IDs which can be used as seeds
    %   for the edge queries
    %
    % skipEdges
    %   Optional Nx2 matrix. Each row represents an edge in
    %   global segment IDs that will be ignored from the edge
    %   sampling procedure.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    batchSize = 250;
    
    % get all edges
    edges = graph.edges;
    
    if exist('skipEdges', 'var')
        % remove edges to skip
        keepMask = ~ismember( ...
            edges, skipEdges, 'rows');
        edges = edges(keepMask, :);
    end
    
    % find edges to / from 'segIds'
    keepMask = any(ismember(edges, segIds), 2);
    edges = edges(keepMask, :);
    
    % pick a random batch
    edgeCount = size(edges, 1);
    pickCount = min(batchSize, edgeCount);
    
    % prepare output
    pickIds = randperm(edgeCount, pickCount);
    edges = edges(pickIds, :);
end
