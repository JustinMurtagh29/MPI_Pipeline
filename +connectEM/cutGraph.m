function graphCut = cutGraph(p, graph, segmentMeta, borderMeta, heuristics, ...
        borderSizeThreshold, segmentSizeThreshold)
    % Restrict graph based on heuristics results and border and segment size threshold

    % Hack number 42 part 10: Exclude all ER components as well
    er = connectEM.getERcomponents();
    heuristics.mapping = cat(2, heuristics.mapping, {cat(1, er{:})});
    heuristicIdx(cat(1, er{:})) = true;

    % Keep only edges above borderSizeThreshold (and correspondences)
    corrIdx = isnan(graph.borderIdx);
    edgeIdx = false(size(corrIdx));
    edgeIdx(corrIdx) = true;
    borderSizes = borderMeta.borderSize(graph.borderIdx(~corrIdx));
    edgeIdx(~corrIdx) =  borderSizes > borderSizeThreshold;
    % Hack 42 part 11:
    edgeIdx(corrIdx) = false;
    remainingEdges = graph.edges(edgeIdx, :);
    remainingProb = graph.prob(edgeIdx);
    % Calculate maximum probability remaining for each segment and exclude based on both thresholds
    maxProb = accumarray(cat(1,remainingEdges(:,1),remainingEdges(:,2)), cat(1,remainingProb, remainingProb),[segmentMeta.maxSegId 1], @max);
    smallIdx = segmentMeta.voxelCount <= segmentSizeThreshold & ~heuristics.heuristicIdx;
    lowProbIdx = segmentMeta.voxelCount > segmentSizeThreshold & maxProb <= 0.5 & ~heuristics.heuristicIdx;

    % Remove heuristics and small or 'disconnected' segments from graph
    removedIds = cat(1, heuristics.mapping{:}, find(smallIdx), find(lowProbIdx));
    keptIds = setdiff(1:double(segmentMeta.maxSegId), removedIds);
    keepEdgeIdx = all(ismember(remainingEdges, keptIds), 2);
    graphCut.edges = remainingEdges(keepEdgeIdx,:);
    graphCut.prob = remainingProb(keepEdgeIdx);

end

