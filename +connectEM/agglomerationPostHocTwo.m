function agglomerationPostHocTwo(options, filename, graph)
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
    if ~exist('graph', 'var')
        graph = load([p.saveFolder 'graph.mat'], 'prob', 'edges', 'borderIdx');
    end
    load([p.saveFolder 'directions.mat']);
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    gridAgglo_05{564} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat');
    segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId', 'cubeIdx');
    segmentMeta.point = segmentMeta.point';
    segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

    % find relevant edges in graph.edge
    borderIdxs = graph.borderIdx(directions.edgeposition);
    borderSizeFake = inf(size(directions.edgeposition));
    borderSizeFake(~isnan(borderIdxs)) = borderMeta.borderSize(borderIdxs(~isnan(borderIdxs)));
    forceKeepEdges = borderSizeFake > options.borderSizeThreshold;
    forceKeepEdges = forceKeepEdges & graph.prob(directions.edgeposition) > options.probThreshold;
    forceKeepEdges = forceKeepEdges & directions.scores > options.scoreThreshold;
    forceKeepEdges = forceKeepEdges & directions.latent(graph.edges(directions.edgeposition, 1))' > options.latentThreshold;
    forceKeepEdges = forceKeepEdges & segmentMeta.voxelCount(graph.edges(directions.edgeposition, 2)) > options.sizeThreshold;
    forceKeepEdges = forceKeepEdges & segmentMeta.axonProb(graph.edges(directions.edgeposition, 2)) > options.axonProbThreshold;
    forceKeepEdges = forceKeepEdges & directions.agglomerationSize(graph.edges(directions.edgeposition, 1))' > options.agglomerationSizeThreshold;

    % force graph probability
    graph.prob(directions.edgeposition(forceKeepEdges)) = 17;
    % force axon probability
    segmentMeta.axonProb(directions.edges(forceKeepEdges, :)) = 18;
    gridAgglo_05{564}.axonProbThreshold = options.axonProbThreshold;

    optional.forceKeepEdges = directions.edgeposition(forceKeepEdges);
    optional.segmentMeta = segmentMeta;
    optional.borderMeta = borderMeta;
    optional.skipDendrites = true;
    optional.calculateMetrics = true;
    connectEM.agglomerationModify(gridAgglo_05{564}, filename, graph, optional);
end
function donothing()
end
