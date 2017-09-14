function [ids, cEdges] = agglomerateSeed( seedIds, edges, borderSize, borderIdx, mergeP, tp, ts)
    %AGGLOMERATESEED Agglomerate from a seed segment.
    % INPUT seedIds: [Nx1] int
    %           Segmentation ids of the starting segments.
    %       edges: [Nx2] int
    %           Edge list.
    %       borderSize: [Nx1] double
    %           BorderSize list.
    %       borderIdx: [Nx1] double
    %           BorderIdx list.
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

    cEdges = edges(mergeP > tp, :);
    cBorderIdx = borderIdx(mergeP > tp, :);
    cBorderSize = borderSize(mergeP > tp, :);
    cEdges = cEdges(cBorderSize > ts, :);
    
    disp('findConnectedComponents');
    cc = Graph.findConnectedComponents(cEdges, false, false);
    idx = find(cellfun(@(x)any(ismember(seedIds, x)), cc));
    ids = cell(1,1);
    for i=1:size(idx,1)    
        if size(cc{idx(i)},1) > 25
            ids{i} = cc{idx(i)}';
        end
    end
    ids = cell2mat(ids)';
    ids = unique(ids);

end
