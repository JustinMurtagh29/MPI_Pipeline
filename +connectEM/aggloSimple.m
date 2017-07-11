function aggloSimple(borderSize, probThreshold, sizeThreshold, outputFolder, optional);

    borderSize = 100;
    probThreshold = .99;
    sizeThreshold = 500;

    % Start by loading parameter file
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    % To keep workspace clean here remove parameter for training (pT)
    clear pT;
    if ~exist('optional', 'var')
        optional = [];
    end

    display('Loading data:');
    tic;
    % Load new graph representation including (resorted) borderIdx [p.saveFolder 'graphNew.mat']
    % These indicate which index in [p.saveFolder 'globalBorder.mat'] each edge corresponds to
    % Correspondences have NaN as borderIdx
    % Load 'neighbours' and 'neighProb' in addition if you want to do (many) local searches in the graph
    if ~exist('graph', 'var')
        graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges', 'borderIdx');
    end
    % Load information about edges
    if isfield(optional, 'borderMeta')
        borderMeta = optional.borderMeta;
    else
        borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
    end
    % Load meta information of segments
    if isfield(optional, 'segmentMeta')
        segmentMeta = optional.segmentMeta;
    else
        segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId', 'cubeIdx');
        segmentMeta.point = segmentMeta.point';
    end
    toc;

    display('Loading data:');
    tic;
    % Keep only edges above borderSizeThreshold (and correspondences)
    corrIdx = isnan(graph.borderIdx);
    edgeIdx = false(size(corrIdx));
    edgeIdx(corrIdx) = true;
    borderSizes = borderMeta.borderSize(graph.borderIdx(~corrIdx));
    edgeIdx(~corrIdx) =  borderSizes > borderSizeThreshold;
    graphCut.edges = graph.edges(edgeIdx,:);
    graphCut.prob = graph.prob(edgeIdx);
    toc;

    display('Performing agglomeration on graph:');
    tic;
    [agglos, aggloSize, aggloEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCut, segmentMeta, probThreshold, sizeThreshold);
    toc;

    display('Writing skeletons for debugging the process:');
    tic;
    connectEM.skeletonFromAgglo(graphCut.edges, segmentMeta, ...
        agglos, 'agglos', outputFolder);
    toc;

end

