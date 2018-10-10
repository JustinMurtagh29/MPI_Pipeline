function [mergedBorders, mergeScores] = run(class, edges, borders)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    edgeT = table;
    edgeT.edge = edges;
    edgeT.borderIds = num2cell(transpose(1:size(edgeT, 1)));
    edgeT.score = cellfun(@(ids) median(class( ...
        unique(cell2mat(borders(ids))))), edgeT.borderIds);

    segCount = numel(unique(edges));
    mergedBorders = cell(segCount - 1, 1);
    mergeScores = nan(segCount - 1, 1);

    for curId = 1:(segCount - 1)
       [curScore, curEdgeId] = max(edgeT.score);
        curEdge = edgeT.edge(curEdgeId, :);

        mergedBorders(curId) = edgeT.borderIds(curEdgeId, :);
        mergeScores(curId) = curScore;
        
        curMask = ismember(edgeT.edge, curEdge);
        edgeT.edge(curMask) = curEdge(1);
        edgeT.edge = sort(edgeT.edge, 2);
        edgeT(all(curMask, 2), :) = [];
        if isempty(edgeT); break; end
        
        %% Update interfaces
       [~, curUni, curOldToUni] = unique(edgeT.edge, 'rows');
        curUniCount = accumarray(curOldToUni, 1);

        % NOTE(amotta): Which edges consist of multiple borders? The scores
        % of these edges need to be updated by merging the borders and then
        % recalculating the score.
        curMultUni = curUni(curUniCount > 1);
        if isempty(curMultUni); continue; end
       [~, curOldToUni] = ismember(curOldToUni, find(curUniCount > 1));

        % Merge borders and update scores for edges with multiplicity.
        edgeT.borderIds(curMultUni) = accumarray( ...
            nonzeros(curOldToUni), find(curOldToUni), [], ...
            @(ids) {cell2mat(edgeT.borderIds(ids))});
        edgeT.score(curMultUni) = cellfun(@(ids) median( ...
            class(unique(cell2mat(borders(ids))))), ...
            edgeT.borderIds(curMultUni));
        edgeT = edgeT(curUni, :);
    end
    
    % Sanity check
    assert(isequal(isnan(mergeScores), ...
        cellfun(@isempty, mergedBorders)));
    
    mergedBorders(isnan(mergeScores)) = [];
    mergeScores(isnan(mergeScores)) = [];
end
