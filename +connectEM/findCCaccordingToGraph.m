function [cc, sizeCC] = findCCaccordingToGraph(graph, ids, segmentMeta)

    % All possible edges between ids
    idx = all(ismember(graph.edges, ids),2);
    edges = graph.edges(idx,:);

    % Find connected components
    cc = Graph.findConnectedComponents(edges, false, true);

    % Keep all that are not connected as a seperate component
    missingIds = setdiff(ids, cat(1,cc{:}));
    cc(end+1:end+length(missingIds),1) = mat2cell(missingIds, ones(length(missingIds),1));

    % Sort according to size of cc
    sizeCC = cellfun(@(x)sum(segmentMeta.voxelCount(x)), cc);
    [sizeCC, idxResort] = sort(sizeCC, 'descend');
    cc = cc(idxResort);

    % Restrict to 1 million voxel for now as well
    idx = sizeCC > 1e6;
    cc = cc(idx);
    sizeCC = sizeCC(idx);
    
end

