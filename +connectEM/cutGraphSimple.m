function [graphCut, proxyFilter] = cutGraphSimple(p, graph, segmentMeta, borderMeta, ...
        borderSizeThreshold, segmentSizeThreshold)
    % Restrict graph based on heuristics results and border and segment size threshold
    % Keep only edges above borderSizeThreshold (and correspondences)
    corrIdx = isnan(graph.borderIdx);
    edgeIdx = false(size(corrIdx));
    edgeIdx(corrIdx) = true;
    borderSizes = borderMeta.borderSize(graph.borderIdx(~corrIdx));
    edgeIdx(~corrIdx) =  borderSizes > borderSizeThreshold;
    forceKeepEdges = corrIdx;% SL temp hack
    %edgeIdx(forceKeepEdges) = true;
    remainingEdges = graph.edges(edgeIdx, :);
    remainingProb = graph.prob(edgeIdx);
    proxyFilter = false(size(graph.prob));
    proxyFilter(forceKeepEdges) = true;
    proxyFilter = proxyFilter(edgeIdx);

    % Calculate maximum probability remaining for each segment and exclude based on both thresholds
    maxProb = accumarray(cat(1,remainingEdges(:,1),remainingEdges(:,2)), cat(1,remainingProb, remainingProb),[segmentMeta.maxSegId 1], @max);
    smallIdx = segmentMeta.voxelCount <= segmentSizeThreshold;
    lowProbIdx = segmentMeta.voxelCount > segmentSizeThreshold & maxProb <= 0.5;

    % Remove heuristics and small or 'disconnected' segments from graph (and cubes excluded due to catastrpohic merger/alignment)
    removedIds = cat(1, find(smallIdx), find(lowProbIdx));
    keptIds = setdiff(1:double(segmentMeta.maxSegId), removedIds);
    keepEdgeIdx1 = all(ismember(remainingEdges, keptIds), 2);

    %removedIds = cat(1, heuristics.mapping{:}, find(excludedSegmentIdx));
    %keptIds = setdiff(1:double(segmentMeta.maxSegId), removedIds);
    %keepEdgeIdx2 = all(ismember(remainingEdges, keptIds), 2);
    graphCut.edges = remainingEdges(keepEdgeIdx1 | proxyFilter,:);
    graphCut.prob = remainingProb(keepEdgeIdx1 | proxyFilter);

end
