function gtConsolidated = consolidateRedundantAnnotations( p, gt )
% Consolidate redundant annotations as extracted by
% getContinuityLabelsFromNml

gtConsolidated = struct();
for i=1:size(gt,1)
    allEdges = cat(1,gt(i,:).edges);
    allLabels = cat(1,gt(i,:).labels);
    uniqueSegments = unique(allEdges(:));
    [~, allEdgesIdx] = ismember(allEdges, uniqueSegments);
    edgeVoteMatrix = accumarray(allEdgesIdx, allLabels, repmat(length(uniqueSegments),1,2));
    % Need to add all correspondences here, otherwise splits up skeletons
    corr = Seg.Global.getGlobalCorrespondences(p);
    idx = all(ismember(corr, uniqueSegments),2);
    [~, allCorrIdx] = ismember(corr, uniqueSegments);
    corrMatrix = accumarray(allCorrIdx(idx,:), 4, repmat(length(uniqueSegments),1,2));
    edgeVoteMatrix = edgeVoteMatrix + corrMatrix;
    [edges(:,1), edges(:,2)] = find(edgeVoteMatrix >= 2);
    edges = changem(edges, uniqueSegments, 1:length(uniqueSegments));
    gtConsolidated(i,1).edges = edges;
    gtConsolidated(i,1).segIds = Graph.findConnectedComponents(edges, false, true);
    clear edges;
end

end

% This function is supposed to be present in Matlab Mapping Toolbox, just
% reimplemented here as I was not able to checkout the toolbox
function newmap = changem(map, newcode, oldcode)

    [~, idx] = ismember(map, oldcode);
    newmap = newcode(idx);

end
