
    % Parameters
    outputFolder = '/gaba/scratch/mberning/axonQueryGeneration/';
    minAxonLength = 5000; % use only these for query generation and target in remapped segmentation
    bboxDist = 1000;

    % Load data
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    [graph, segmentMeta, borderMeta, globalSegmentPCA] = connectEM.loadAllSegmentationData(p);

    % Load state of axon agglomeration
    load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat', 'axonsNew');

    % Calculate axon directionality
    options.voxelSize = p.raw.voxelSize;
    options.bboxDist = bboxDist;
    directionality = connectEM.calculateDirectionalityOfAgglomerates(axonsNew, graph, segmentMeta, borderMeta, globalSegmentPCA, options);
    save([outputFolder 'directionality.mat'], 'directionality', '-v7.3');

    % Calculate size as maximal distance between all agglomerates in an agglomeration
    calculateLength = @(x)max(pdist(bsxfun(@times, double(borderMeta.borderCoM(x, :)), p.raw.voxelSize)));
    axonLength = cellfun(calculateLength, directionality.borderIdx, 'uni', 0);
    save([outputFolder 'axonLength.mat'], 'axonLength');

    % Write out skeletons for current agglo
    segmentMeta.point = segmentMeta.point';
    y = connectEM.evaluateAggloMetaMeta(graph, axonsNew, [], 'agglosForAxonQueryGeneration', segmentMeta);
    segmentMeta.point = segmentMeta.point';

    % Keep only those axons and directionality metrics that are longer that specified minimal length
    idxEmpty = cellfun('isempty', axonLength);
    idxKeep = ~idxEmpty;
    idxKeepLarge = cellfun(@(x)x > minAxonLength, axonLength(idxKeep));
    idxKeep(idxKeep) = idxKeepLarge;
    % Save small agglomerates seperatly for possible use later
    axonsSmall = axonsNew(~idxKeep);
    save([outputFolder 'axonsSmall.mat'], 'axonsSmall');
    axonsNew = axonsNew(idxKeep);
    directionality = structfun(@(x)x(idxKeep), directionality, 'uni', 0);
    axonLength = cell2mat(axonLength(idxKeep));
    clear idx*;

    % Write new segmentation based on axon queries
    mapping = connectEM.createLookup(segmentMeta, axonsNew);
    Seg.Global.applyMappingToSegmentation(p, mapping, [outputFolder '1/']);
    clear mapping;

    % Save state for bookkeeping
    save([outputFolder 'beforeQueryGeneration.mat'], 'directionality', 'axonsNew', 'axonLength', '-v7.3');

    % Generate axon queries
    querySaveFolder = [outputFolder 'queriesMat/'];
    if ~exist(querySaveFolder, 'dir')
        mkdir(querySaveFolder)
    end
    connectEM.generateAxonQueries(p, graph, segmentMeta, borderMeta, directionality, axonsNew, querySaveFolder);

