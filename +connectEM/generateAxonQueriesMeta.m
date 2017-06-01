
    % Parameters
    outputFolder = '/gaba/scratch/mberning/axonQueryGeneration/';
    minAxonLength = 5000; % use only these for query generation and target in remapped segmentation
    bboxDist = 1000;

    % Load data
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    [graph, segmentMeta, borderMeta, globalSegmentPCA] = loadAllSegmentationData(p);

    % TODO: Decide based on which agglo to use once grid search is finished
    % Load state of axon agglomeration
    load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00025/metricsFinal.mat', 'axonsNew');

    % Calculate axon directionality
    options.voxelSize = p.raw.voxelSize;
    options.bboxDist = bboxDist; 
    directionality = connectEM.calculateDirectionalityOfAgglomerates(axonsNew, graph, segmentMeta, borderMeta, globalSegmentPCA, options);
    save([outputFolder 'directionality.mat'], 'directionality');

    % Calculate size as maximal distance between all agglomerates in an agglomeration
    calculateLength = @(x)max(pdist(bsxfun(@times, double(borderMeta.borderCoM(x, :)), p.raw.voxelSize)));
    axonLength = calculateLength(directionality.borderIdx);
    
    % Keep only those axons and directionality metrics that are longer that specified minimal length
    idx = axonLength > minAxonLength;
    axonsNew = axonsNew(idx);
    directionality = cellfun(@(x)x(idx), directionality, 'uni', 0);
    clear idx;

    % Write new segmentation based on axon queries
    mapping = connectEM.createLookup(segmentMeta, axonsNew);
    Seg.Global.applyMappingToSegmentation(p, mapping, [outputFolder '1/']);

    % Generate axon queries 
    connectEM.generateAxonQueries(p, graph, segmentMeta, borderMeta, directionality, axonsNew);

    % Visualize axon queries as skeletons for debugging process
    idx = randperm(numel(axonsNew), 100);
    connectEM.debugNewQueries(segmentMeta, axonsNew(idx), q(idx), [outputFolder 'queryVisualization/']);
    clear idx;

