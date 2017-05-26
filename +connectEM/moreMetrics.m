function metrics = moreMetrics(axonsNew, name, segmentMeta)
    
    voxelCount = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axonsNew);
    [metrics.voxelCount, idx] = sort(voxelCount, 'descend');
    axonsNew = axonsNew(idx);
    clear voxelCount idx;
    % Calculate metrics as before
    metrics.y = connectEM.evaluateAggloMetaMeta(graph, axonsNew, [], name, segmentMeta);
    % And only on agglomerates with more than 5 micron path length (if no large axons percolators, those runs are bad anyways)
    if metrics.y.axonPercolators(1) > 1e8
        metrics.percolation = true;
    else
        metrics.percolation = false;
        metrics.pathLength  = connectEM.getPathLengthFromAgglomeration(axonsNew, segmentMeta.point);
        metrics.nrAgglos = numel(metrics.pathLength);
        metrics.pathLengthAgglos = sum(metrics.pathLength);
        metrics.nrLargeAgglos = sum(metrics.pathLength > 5);
        metrics.pathLengthLargeAgglos = sum(metrics.pathLength(metrics.pathLength > 5));
        metrics.yLarge = connectEM.evaluateAggloMetaMeta(graph, axonsNew(metrics.pathLength > 5), [], name, segmentMeta);
    end

end

