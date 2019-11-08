function [partition, partitionSize] = partitionSortAndKeepOnlyLarge(graph, segmentMeta, probThreshold, sizeThreshold, corrEdges, maxSegId)
    % Partition graph based on probThreshold, sort CC by voxel size and keep only larger than sizeThreshold in voxel

    [partition, partitionEdges] = connectEM.partitionWholeDataset(graph, probThreshold);
    
    % force agglos that are at correspondances to be combined
    [~,partition] = Agglo.mergeAgglosAtEdges(partition, corrEdges, maxSegId);

    % Sort agglomerates by voxel size
    partitionSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), partition);
    [partitionSize, idx] = sort(partitionSize, 'descend');
    partition = partition(idx);
    % Keep only agglomerates that have at least sizeThreshold million voxels
    idx = partitionSize > sizeThreshold;
    partition = partition(idx);
    partitionSize = partitionSize(idx);
    clear idx;

end

