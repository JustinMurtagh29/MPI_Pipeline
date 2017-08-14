function collectGlobalGraphStruct(p)
    edgeFile = fullfile(p.saveFolder, 'globalEdges.mat');
    probFile = fullfile(p.saveFolder, 'globalNeuriteContinuityProb.mat');
    corrFile = fullfile(p.saveFolder, 'mapping.mat');

    % The edge list only contains edges between segments in the same
    % segmentation cube. Edges between different cubes are called
    % 'correspondences' and must be added separately here.
    edges = load(edgeFile, 'edges');
    prob = load(probFile, 'prob');
    
    corr = load(corrFile, 'correspondences');
    corr = corr.correspondences;
    
    borderIdx = [
        reshape(1:size(edges, 1), [], 1); % for border-based edges
        nan(size(corr, 1), 1)];           % for correspondences
    
    edges = [edges; corr];
    prob = [prob; ones(size(corr, 1), 1)];

    % sort edges and probabilities
   [edges, rows] = sortrows(edges);
    prob = prob(rows);
    borderIdx = borderIdx(rows);

    % build output
    graph = struct;
    graph.edges = edges;
    graph.prob = prob;
    graph.borderIdx = borderIdx;

    % save output
    graphFile = fullfile(p.saveFolder, 'graph.mat');
    Util.saveStruct(graphFile, graph);
end

