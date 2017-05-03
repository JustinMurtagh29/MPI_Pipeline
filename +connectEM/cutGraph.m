function graphCut = cutGraph(p, graph, segmentMeta, borderMeta, heuristics, ...
        borderSizeThreshold, segmentSizeThreshold, forceKeepEdges)
    % Restrict graph based on heuristics results and border and segment size threshold

    % Still hacky, but keep for now: Exclude all ER components and all cube containing catastrphic merger
    [er, cm] = connectEM.getERcomponents();
    heuristics.mapping = cat(2, heuristics.mapping, {cat(1, er{:})});
    heuristicIdx(cat(1, er{:})) = true;
    excludedCubeIdx = unique(cellfun(@(x)mode(segmentMeta.cubeIdx(x)), cm));
    excludedSegmentIdx = ismember(segmentMeta.cubeIdx, excludedCubeIdx);

    % Keep only edges above borderSizeThreshold (and correspondences)
    corrIdx = isnan(graph.borderIdx);
    edgeIdx = false(size(corrIdx));
    edgeIdx(corrIdx) = true;
    borderSizes = borderMeta.borderSize(graph.borderIdx(~corrIdx));
    edgeIdx(~corrIdx) =  borderSizes > borderSizeThreshold;
    edgeIdx(forceKeepEdges) = true;
    remainingEdges = graph.edges(edgeIdx, :);
    remainingProb = graph.prob(edgeIdx);
    % Calculate maximum probability remaining for each segment and exclude based on both thresholds
    maxProb = accumarray(cat(1,remainingEdges(:,1),remainingEdges(:,2)), cat(1,remainingProb, remainingProb),[segmentMeta.maxSegId 1], @max);
    smallIdx = segmentMeta.voxelCount <= segmentSizeThreshold & ~heuristics.heuristicIdx;
    lowProbIdx = segmentMeta.voxelCount > segmentSizeThreshold & maxProb <= 0.5 & ~heuristics.heuristicIdx;

    % Remove heuristics and small or 'disconnected' segments from graph (and cubes excluded due to catastrpohic merger/alignment)
    removedIds = cat(1, heuristics.mapping{:}, find(smallIdx), find(lowProbIdx), find(excludedSegmentIdx));
    keptIds = setdiff(1:double(segmentMeta.maxSegId), removedIds);
    keepEdgeIdx = all(ismember(remainingEdges, keptIds), 2);
    keepEdgeIdx(forceKeepEdges) = true;
    graphCut.edges = remainingEdges(keepEdgeIdx,:);
    graphCut.prob = remainingProb(keepEdgeIdx);

end
