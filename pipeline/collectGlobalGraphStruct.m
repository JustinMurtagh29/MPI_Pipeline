function collectGlobalGraphStruct(p)
    load([p.saveFolder 'globalEdges.mat']);
    load([p.saveFolder 'globalGPProbList.mat']);
    borderIdx = (1:numel(prob))';

    % The edge list only contains edges between segments in
    % the same segmentation cube. Edges between different
    % cubes are called 'correspondences' and must be added
    % separately here.
    % Now using the "new" correspondences, commented out old approach
     corrEdges = Seg.Global.getGlobalCorrespondences(p);
    % See pipeline/correspondenceTests.m for details on calculation or:
    % https://gitlab.mpcdf.mpg.de/connectomics/pipeline/commit/4b20a81f7daa29dc6fa555913afb265f538fafe9
    % Loads "corrEdges"

    corrProb  = ones(size(corrEdges, 1), 1);
    corrBorderIdx = NaN(size(corrEdges, 1), 1);

    % build total edges and prob fields
    edges = [edges; corrEdges];
    prob  = [prob ; corrProb ];
    borderIdx = [borderIdx; corrBorderIdx];

    % sort edges and probabilities
    [edges, edgeRows] = sortrows(edges);
    prob = prob(edgeRows);
    borderIdx = borderIdx(edgeRows);

    % build look up tables for neighbouring
    % segment IDs and their connection probabilities
    [neighbours, neighboursIdx] = ...
        Graph.edges2Neighbors(edges);
    neighProb = cellfun(@(x) ...
        prob(x), neighboursIdx, 'UniformOutput', false);
    neighBorderIdx = cellfun(@(x) ...
        borderIdx(x), neighboursIdx, 'UniformOutput', false);

    % build output
    graph = struct;
    graph.neighbours = neighbours;
    graph.neighProb = neighProb;
    graph.neighBorderIdx = neighBorderIdx;
    graph.edges = edges;
    graph.prob = prob;
    graph.borderIdx = borderIdx;

    % save output
    graphFile = [p.saveFolder, 'graphNew.mat'];
    Util.saveStruct(graphFile, graph);
end

