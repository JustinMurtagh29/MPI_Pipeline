function forcingNum = agglomerationPostHocTwo(options, filename, graph, borderMeta, segmentMeta, directions, agglo, idx,topfolder)
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
    forceCorrespondences = isnan(graph.borderIdx) & any(segmentMeta.axonProb(graph.edges) > 0.5, 2);

    % force graph probability
    oldForce = load([topfolder, 'forceKeepEdges_' num2str(idx - 1, '%0.3u')], 'forceKeepEdgesStore')
    forceKeepEdgesStore = [oldForce.forceKeepEdgesStore; directions.edgeposition(forceKeepEdges)];
    save([topfolder, 'forceKeepEdges_' num2str(idx, '%0.3u')], 'forceKeepEdgesStore');

    graph.prob(forceKeepEdgesStore) = 17;
    graph.prob(forceCorrespondences) = 19;
    % force axon probability
    segmentMeta.axonProb(graph.edges(forceKeepEdgesStore, :)) = 18;
    segmentMeta.axonProb(graph.edges(forceCorrespondences, :)) = 20;
    agglo.axonProbThreshold = options.axonProbThreshold;
    optional.forceKeepEdges = [forceKeepEdgesStore; find(forceCorrespondences)];
    optional.segmentMeta = segmentMeta;
    optional.borderMeta = borderMeta;
    optional.skipDendrites = true;
    optional.calculateMetrics = false;

    connectEM.agglomerationModify(agglo, filename, graph, optional);
    forcingNum = sum(forceKeepEdges);
end
function donothing()
end
