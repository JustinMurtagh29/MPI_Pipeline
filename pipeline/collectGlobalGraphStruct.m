function collectGlobalGraphStruct (p)

    % The edge list only contains edges between segments in
    % the same segmentation cube. Edges between different
    % cubes are called 'correspondences' and must be added
    % separately here.
    corrEdges = Seg.Global.getGlobalCorrespondences(p);
    corrProb  = ones(size(corrEdges, 1), 1);

    % build total edges and prob fields
    m = load([p.saveFolder 'globalEdges.mat']);
    edges = [m.edges; corrEdges];
    m = load([p.saveFolder 'globalGPProbList.mat']);
    prob  = [m.prob ; corrProb ];
    clear m
    
    % sort edges and probabilities
    [edges, edgeRows] = sortrows(edges);
    prob = prob(edgeRows);

    % build look up tables for neighbouring
    % segment IDs and their connection probabilities
    [neighbours, neighboursIdx] = ...
        Graph.edges2Neighbors(edges);
    neighProb = cellfun(@(x) ...
        prob(x), neighboursIdx, 'UniformOutput', false);

    % build output
    graph = struct;
    graph.neighbours = neighbours;
    graph.neighProb = neighProb;
    graph.edges = edges;
    graph.prob = prob;

    % save output
    graphFile = [p.saveFolder, 'graph.mat'];
    Util.saveStruct(graphFile, graph);
end

