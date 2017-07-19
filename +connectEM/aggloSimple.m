function aggloSimple(p,borderSizeThreshold, probThreshold, sizeThreshold, outputFolder, optional)
    if nargin < 2 || isempty(borderSizeThreshold)
       borderSizeThreshold = 100;
    end
    if nargin < 3 || isempty(probThreshold)
       probThreshold = .995;
    end
    if nargin < 4 || isempty(sizeThreshold)
       sizeThreshold = 1000;
    end
    if nargin < 1 || isempty(p)
       % Start by loading parameter file
       load('/gaba/scratch/mbeining/testseg_KSMBtmp/PC4/allParameter.mat');
    end
    % To keep workspace clean here remove parameter for training (pT)
    clear pT;
    if ~exist('optional', 'var')
        optional = [];
    end


    parameters.scale.x = p.raw.voxelSize(1);
    parameters.scale.y = p.raw.voxelSize(2);
    parameters.scale.z = p.raw.voxelSize(3);
    parameters.offset.x = '0';
    parameters.offset.y = '0';
    parameters.offset.z = '0';
    parameters.experiment.name = p.experimentName;




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

    display('Cutting graph:');
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
    % Use only agglos 2:end here as first is large percolator, you should try to find out what/why that is
    connectEM.skeletonFromAgglo(graphCut.edges, segmentMeta, ...
        agglos(2:end), 'agglos', outputFolder,parameters);
    toc;

end

