% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% Load membrane and segmentation data
load('/home/amotta/Desktop/mrnet.mat');

%% Build edges and borders
[edges, borders] = SynEM.Svg.findEdgesAndBorders(seg);
borders = reshape({borders.PixelIdxList}, [], 1);

%%
edgeT = table;
edgeT.edge = edges;
edgeT.vxIds = borders;
edgeT.score = cellfun(@(ids) median(class(ids)), edgeT.vxIds);

maxSegId = max(edgeT.edge(:));
mergedEdges = nan(maxSegId - 1, 3);
for curId = 1:(maxSegId - 1)
   [~, curUni, curOldToUni] = unique(edgeT.edge, 'rows');
    curUniCount = accumarray(curOldToUni, 1);
    
    % NOTE(amotta): Which edges consist of multiple borders? The scores of
    % these edges need to be updated by merging the borders and then
    % recalculating the score.
    curMultUni = curUni(curUniCount > 1);
   [~, curOldToUni] = ismember(curOldToUni, find(curUniCount > 1));
    
    % Merge borders and update scores for edges with multiplicity.
    edgeT.vxIds(curMultUni) = accumarray( ...
        nonzeros(curOldToUni), find(curOldToUni), [], ...
        @(ids) {cell2mat(edgeT.vxIds(ids))});
    edgeT.score(curMultUni) = cellfun( ...
        @(ids) median(class(ids)), ...
        edgeT.vxIds(curMultUni));
    edgeT = edgeT(curUni, :);
    
   [curMaxScore, curEdgeId] = max(edgeT.score);
    curEdge = edgeT.edge(curEdgeId, :);
    
    mergedEdges(curId, 1:2) = curEdge;
    mergedEdges(curId, 3) = curMaxScore;
    
    % Assume to merge
    curMask = ismember(edgeT.edge, curEdge);
    
    edgeT.edge(curMask) = maxSegId + 1;
    edgeT.edge = sort(edgeT.edge, 2);
    edgeT(all(curMask, 2), :) = [];
    
    maxSegId = maxSegId + 1;
end
