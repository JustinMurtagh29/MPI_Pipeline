function [graphCut, mapping, mappingNames, mappingSize] = cutGraph(p, graph, segmentMeta, borderMeta, ...
        borderSizeThreshold, segmentSizeThreshold)
    % Restrict graph based on heuristics results

    %{
    % Segments classified by heuristics
    load([p.saveFolder 'heuristicResult.mat']);
    assert(length(segIds) == max(segIds));
    % Vessel and endothelial cells
    vesselIdx = vesselScore > 0.5;
    %endoIdx = growOutHeuristics(graph, segmentMeta, vesselIdx, 0.995, 1000);
    temp = load([p.saveFolder 'perivessel.mat']);
    endoIdx = false(segmentMeta.maxSegId, 1);
    endoIdx(cat(1, temp.comps{:})) = true;
    endoIdx = endoIdx & ~vesselIdx;
    clear temp;
    % Nuclei + added if not grown to completion
    nucleiIdx = nucleiScore > 0.5 & ~vesselIdx & ~endoIdx;
    addedNucleiIdx = growOutHeuristics(graph, segmentMeta, nucleiIdx, 0.999, 1000);
    addedNucleiIdx = addedNucleiIdx & ~ vesselIdx & ~endoIdx;
    % Myelin + added to try to avoid merger with myelin
    myelinIdx = myelinScore > 0.5 & ~vesselIdx & ~endoIdx & ~nucleiIdx & ~addedNucleiIdx;
    addedMyelinIdx = growOutHeuristics(graph, segmentMeta, myelinIdx, 0.995, 1000);
    addedMyelinIdx = addedMyelinIdx & ~vesselIdx & ~endoIdx & ~nucleiIdx & ~addedNucleiIdx;
    % All heuristics exclusions
    heuristicIdx = vesselIdx | endoIdx | nucleiIdx | addedNucleiIdx | myelinIdx | addedMyelinIdx;

    % Keep only segments larger than segmentSizeThreshold voxel and ...
    % that have more than 50% connection probability to another segment after removing border smaller than borderSizeThreshold
    load([p.saveFolder  'globalGPProbList.mat']);
    load([p.saveFolder 'globalEdges.mat']);
    corrEdges = Seg.Global.getGlobalCorrespondences(p);
    corrProb  = ones(size(corrEdges, 1), 1);
    %}
    % Load saved state as this will not change
    load([p.saveFolder 'graphCut.mat']);

    edgeIdx = find(borderMeta.borderSize > borderSizeThreshold);
    remainingEdges = edges(edgeIdx, :);
    remainingProb = prob(edgeIdx);
    % Add correspondences
    remainingEdges = [remainingEdges; corrEdges];
    remainingProb  = [remainingProb; corrProb];
    % Calculate maximum probability remaining for each segment and exclude based on both thresholds
    maxProb = accumarray(cat(1,remainingEdges(:,1),remainingEdges(:,2)), cat(1,remainingProb, remainingProb),[segmentMeta.maxSegId 1], @max);
    smallIdx = segmentMeta.voxelCount <= segmentSizeThreshold & ~heuristicIdx;
    lowProbIdx = segmentMeta.voxelCount > segmentSizeThreshold & maxProb <= 0.5 & ~heuristicIdx;

    % Concatenate exclusions for output
    mapping = {vesselIdx endoIdx nucleiIdx addedNucleiIdx myelinIdx addedMyelinIdx smallIdx lowProbIdx};
    mapping = cellfun(@find, mapping, 'uni', 0);
    mappingNames = {'vessel' 'endothelial' 'nuclei' 'addedNuclei' 'myelin' 'addedMyelin' 'smallSegments' 'disconnectedSegments'};
    mappingSize = cellfun(@numel, mapping);

    % Remove heuristics and small or disconnected segments from graph
    removedIds = cat(1, mapping{:});
    keptIds = setdiff(1:double(segmentMeta.maxSegId), removedIds);
    keepEdgeIdx = all(ismember(remainingEdges, keptIds), 2);
    graphCut.edges = remainingEdges(keepEdgeIdx,:);
    graphCut.prob = remainingProb(keepEdgeIdx);
    % Sort edges again (e.g. make correspondences find their place)
    [graphCut.edges, edgeRows] = sortrows(graphCut.edges);
    graphCut.prob = graphCut.prob(edgeRows);

end

function addedSegmentIdx = growOutHeuristics(graph, segmentMeta, initialSegmentIdx, probThreshold, sizeThreshold)
    % Grow out heuristics to avoid merger into endothelial, missed myelin etc.

    partition = connectEM.partitionSortAndKeepOnlyLarge(graph, segmentMeta, probThreshold, sizeThreshold);
    initialSegmentIds = find(initialSegmentIdx);
    sizeComponent = cellfun(@numel, partition);
    initialSegmentsPerComponent = cellfun(@(x)sum(ismember(x, initialSegmentIds)), partition);
    initialSegmentFraction = initialSegmentsPerComponent ./ sizeComponent;
    allSegmentIdx = false(segmentMeta.maxSegId, 1);
    allSegmentIdx(cat(1, partition{initialSegmentFraction > 0.1})) = true;
    addedSegmentIdx = allSegmentIdx & ~initialSegmentIdx;

end

