function [ids, cEdges] = agglomerateSeed( seedIds, edges, borderSize, ...
    mergeP, tp, ts)
    %AGGLOMERATESEED Agglomerate from a seed segment.
    % INPUT seedIds: [Nx1] int
    %           Segmentation ids of the starting segments.
    %       edges: [Nx2] int
    %           Edge list.
    %       borderSize: [Nx1] double
    %           BorderSize list.
    %       mergeP: [Nx1] double
    %           Merge probability for the corresponding edges.
    %       tp: double
    %           Lower probability threshold for merging.
    %       ts: double
    %           Lower size threshold for merging.
    % OUTPUT ids: [Nx1] int
    %           Ids that are connected to the seedId.
    %        cEdges: [Nx2] int
    %           The edges that connect the ids.
    % Author: Benedikt Staffler, modified by Robin Hesse

    toKeepEdges = mergeP > tp;
    cEdges = edges(toKeepEdges, :);
    cBorderSize = borderSize(toKeepEdges, :);
    cEdges = cEdges(cBorderSize > ts, :);

    isSeed = false(max(cEdges(:)), 1);
    isSeed(seedIds) = true;
    cc = Graph.findConnectedComponents(cEdges, false, false);
    idx = cellfun(@(x)any(isSeed(x)), cc);
    ids = cell2mat(cc(idx));
end
