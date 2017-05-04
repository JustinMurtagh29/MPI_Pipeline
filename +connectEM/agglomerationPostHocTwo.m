function agglomerationPostHocTwo(options, filename)
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
    graph = load([p.saveFolder 'graph.mat'], 'prob', 'edges', 'borderIdx');
    directions = load([p.saveFolder 'directions.mat'])
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    gridAgglo_05{564} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat');
    segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId', 'cubeIdx');
    segmentMeta.point = segmentMeta.point';
    segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

    % find relevant edges in graph.edge
    forceKeepEdges = borderMeta.borderSize(graph.borderIdx(direction.edgepositions)) > options.borderSizeThreshold;
    forceKeepEdges = forceKeepEdges & graph.prob(direction.edgespositions) > options.probThreshold;
    forceKeepEdges = forceKeepEdges & directions.scores > options.scoreThreshold;
    forceKeepEdges = forceKeepEdges & directions.latent(graph.edges(direction.edgepositions, 1)) > options.latentThreshold;
    forceKeepEdges = forceKeepEdges & segmentMeta.voxelCount(graph.edges(direction.edgeposition, 2)) > options.sizeThreshold;
    forceKeepEdges = forceKeepEdges & segmentMeta.axonProb(graph.edges(direction.edgeposition, 2)) > options.axonProbThreshold;
    forceKeepEdges = forceKeepEdges & directions.agglomerationSize > options.agglomerationSizeThreshold;

    % force graph probability
    graph.prob(direction.edgespositions(forceKeepEdges)) = 17;
    % force axon probability
    segmentMeta.axonProb(direction.edges(fourceKeepEdges, :)) = 18;
    gridAgglo_05{564}.axonProbThreshold = options.axonProbThreshold;

    optional.forceKeepEdges = direction.edgespositions(forceKeepEdges);
    optional.segmentMeta = segmentMeta;
    optional.borderMeta = borderMeta;
    optional.skipDendrites = true;

    connectEM.agglomerationModify(gridAgglo_05{564}, filename, graph, optional);
end
function donothing()
end
