function graphCut = cutGraph(p, graph, segmentMeta, borderMeta, heuristics, ...
        borderSizeThreshold, segmentSizeThreshold)
    % Restrict graph based on heuristics results and border and segment size threshold

    heuristics = load([p.saveFolder 'graphCutNew.mat']);
    heuristics.mapping = {heuristics.vesselIdx heuristics.endoIdx heuristics.nucleiIdx ...
        heuristics.addedNucleiIdx heuristics.myelinIdx heuristics.addedMyelinIdx};
    heuristics.mapping = cellfun(@find, heuristics.mapping, 'uni', 0);

    % Keep only edges above borderSizeThreshold (and correspondences)
    corrIdx = isnan(graph.borderIdx);
    edgeIdx = false(size(corrIdx));
    edgeIdx(~corrIdx) = borderMeta.borderSize(graph.borderIdx(~corrIdx)) > borderSizeThreshold;
    edgeIdx(corrIdx) = true;
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

