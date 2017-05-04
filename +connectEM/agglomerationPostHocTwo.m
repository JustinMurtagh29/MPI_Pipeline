function agglomerationPostHocTwo(options, filename, graph, borderMeta, segmentMeta, directions, agglo)
    % find relevant edges in graph.edge
    borderIdxs = graph.borderIdx(directions.edgeposition);
    borderSizeFake = inf(size(directions.edgeposition));
    borderSizeFake(~isnan(borderIdxs)) = borderMeta.borderSize(borderIdxs(~isnan(borderIdxs)));
    forceKeepEdges = borderSizeFake > options.borderSizeThreshold;
    forceKeepEdges = forceKeepEdges & graph.prob(directions.edgeposition) > options.probThreshold;
    forceKeepEdges = forceKeepEdges & directions.scores > options.scoreThreshold;
    forceKeepEdges = forceKeepEdges & directions.latent(graph.edges(directions.edgeposition, 1))' > options.latentThreshold;
    forceKeepEdges = forceKeepEdges & segmentMeta.voxelCount(graph.edges(directions.edgeposition, 2)) > options.sizeThreshold;
    forceKeepEdges = forceKeepEdges & all(segmentMeta.axonProb(graph.edges(directions.edgeposition, :)) > options.axonProbThreshold, 2);
    forceKeepEdges = forceKeepEdges & directions.agglomerationSize(graph.edges(directions.edgeposition, 1))' > options.agglomerationSizeThreshold;

    % force correspondences
    forceCorrespondences = isnan(graph.borderIdx) & any(segmentMeta.axonProb(graph.edges) > 0.5, 2)

    % force graph probability
    graph.prob(directions.edgeposition(forceKeepEdges)) = 17;
    graph.prob(forceCorrespondences) = 19;
    % force axon probability
    segmentMeta.axonProb(directions.edges(forceKeepEdges, :)) = 18;
    segmentMeta.axonProb(graph.edges(forceCorrespondences)) = 20;
    agglo.axonProbThreshold = options.axonProbThreshold;

    optional.forceKeepEdges = [directions.edgeposition(forceKeepEdges); find(forceCorrespondences)];
    optional.segmentMeta = segmentMeta;
    optional.borderMeta = borderMeta;
    optional.skipDendrites = true;
    optional.calculateMetrics = false;
    connectEM.agglomerationModify(agglo, filename, graph, optional);
end
function donothing()
end
