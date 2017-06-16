function defineHeuristicSegmentExclusions(p)

    % Segments classified by heuristics
    load([p.saveFolder 'heuristicResult.mat']);
    assert(length(segIds) == max(segIds));
    segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
    graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges');
    % Vessel and endothelial cells
    vesselIdx = vesselScore > 0.5;
    % Replaced by heurtistics written by Alessandro
    % endoIdx = growOutHeuristics(graph, segmentMeta, vesselIdx, 0.995, 1000);
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

    % Concatenate exclusions for output
    mapping = {vesselIdx endoIdx nucleiIdx addedNucleiIdx myelinIdx addedMyelinIdx};
    mapping = cellfun(@find, mapping, 'uni', 0);
    mappingNames = {'vessel' 'endothelial' 'nuclei' 'addedNuclei' 'myelin' 'addedMyelin'};
    mappingSize = cellfun(@numel, mapping);

    % Save for speed purposes (growing out nuclei and myelin is expensive part)
    save([p.saveFolder 'heuristicSegmentExclusions.mat']);

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

