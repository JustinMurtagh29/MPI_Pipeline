
    % Load data
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId', 'cubeIdx', 'point');

    % TODO: Decide based on which agglo to use
    % Load state of axon agglomeration
    load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00025/metricsFinal.mat', 'axonsNew');

    % Sort according to size
    axonSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axonsNew);
    [axonSize, idx] = sort(axonSize, 'descend');
    axonsNew = axonsNew(idx);

    % Remove axons that are 
 
