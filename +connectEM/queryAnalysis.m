function [uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = queryAnalysis(segIds, neighbours, nodes, startNode, options)
    % Keep only segmentation IDs that were in bounding box defined in the tracing and return
    % unique segments hit by 27 neighbourhood of all nodes & number of occurences for each.
    % Also return 27 neighbourhood of start node in the same way
    if ~exist('options', 'var')
        options = [];
    end
    if ~isfield(options, 'bboxSize')
        options.bboxSize = 2000;
    end
    % Restrict to bbox (exclude node set in blacked out region in wK)
    extend = round(options.bboxSize ./ [11.24 11.24 28]);
    minPos = startNode - extend;
    maxPos = startNode + extend;
    % Some statistics on all ids of segment ID at node position and all neighboring voxel
    idx1 = all(bsxfun(@ge, nodes, minPos) & bsxfun(@le, nodes, maxPos),2);
    idx2 = all(bsxfun(@eq, nodes, startNode),2);
    allIds = cat(2,segIds(idx1|~idx2),neighbours(idx1|~idx2,:));
    startIds = cat(2, segIds(idx2), neighbours(idx2,:));
    % Check whether some general assumptions hold true
    assert(length(segIds) == size(neighbours,1));
    assert(size(nodes,1) == size(neighbours,1));
    assert(sum(idx2) == 1);
    % How often does each id (except 0) occur at any node execpt start node
    uniqueSegments.segIds = unique(allIds(:));
    uniqueSegments.segIds(uniqueSegments.segIds == 0) = [];
    uniqueSegments.occurences = histc(allIds(:), uniqueSegments.segIds);
    uniqueSegments = sortAccordingToOccurence(uniqueSegments);
    % Same only for start node
    neighboursStartNode.segIds = unique(startIds(:));
    neighboursStartNode.segIds(neighboursStartNode.segIds == 0) = [];
    neighboursStartNode.occurences = histc(startIds(:), neighboursStartNode.segIds);
    neighboursStartNode = sortAccordingToOccurence(neighboursStartNode);
    % Pass some indices as output arguments for visualization (see debugQueryAttachment)
    nodesExcludedIdx = ~idx1;
    startNodeIdx = idx2;
end

function out = sortAccordingToOccurence(in)
    [out.occurences, idxPerm] = sort(in.occurences, 'descend');
    out.segIds = in.segIds(idxPerm);
end
