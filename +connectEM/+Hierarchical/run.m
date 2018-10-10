% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% Load membrane and segmentation data
load('/home/amotta/Desktop/mrnet.mat');

box = [
    4856, 5367; ...
    3289, 3800; ...
    3969, 4224];

%% Build edges and borders
[edges, borders] = SynEM.Svg.findEdgesAndBorders(seg);
borders = reshape({borders.PixelIdxList}, [], 1);
segPoints = Seg.Local.calcSegmentPoint(seg);

%%
edgeT = table;
edgeT.edge = edges;
edgeT.vxIds = borders;
edgeT.score = cellfun(@(ids) median(class(ids)), edgeT.vxIds);

maxSegId = max(edgeT.edge(:));
mergedEdges = nan(maxSegId - 1, 3);

binEdges = linspace(-1.1, +1.1, 101);
fig = figure();
ax = axes(fig);
hold(ax, 'on');

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
        @(ids) {unique(cell2mat(edgeT.vxIds(ids)))});
    edgeT.score(curMultUni) = cellfun( ...
        @(ids) median(class(ids)), ...
        edgeT.vxIds(curMultUni));
    edgeT = edgeT(curUni, :);
    
    if mod(curId, 500) == 1
        histogram( ...
            ax, edgeT.score, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        % keyboard
    end
    
   [curMaxScore, curEdgeId] = max(edgeT.score)
    curEdge = edgeT.edge(curEdgeId, :);
    
    mergedEdges(curId, 1:2) = curEdge;
    mergedEdges(curId, 3) = curMaxScore;
    
    % Assume to merge
    curMask = ismember(edgeT.edge, curEdge);
    
    edgeT.edge(curMask) = curEdge(1);
    edgeT.edge = sort(edgeT.edge, 2);
    edgeT(all(curMask, 2), :) = [];
    
    maxSegId = maxSegId + 1;
end

%% Show result
curThresh = -0.58;
curThreshId = find(mergedEdges(:, 3) < curThresh, 1);

curLUT = 1:maxSegId;

for curId = 1:(curThreshId - 1)
    curLUT(curLUT == mergedEdges(curId, 2)) = mergedEdges(curId, 1);
end

curLUT = uint32([0, curLUT]);
curSegNew = curLUT(1 + seg);

%%
colors = uint8(255 * hsv(double(maxSegId)));
colors = colors(randperm(size(colors, 1)), :);
colors = [0, 0, 0; colors];

curSeg = colors(1 + seg(:), :);
curSeg = reshape(curSeg, [size(seg), 3]);
curSeg = permute(curSeg, [1, 2, 4, 3]);

curSegNew = colors(1 + curSegNew(:), :);
curSegNew = reshape(curSegNew, [size(seg), 3]);
curSegNew = permute(curSegNew, [1, 2, 4, 3]);

%% Build skeleton
skel = skeleton();
skel = skel.setParams( ...
    '2018-10-01_ex144_st08x2', [11.24, 11.24, 28], [0, 0, 0]);

curIds = ceil(linspace(1, size(mergedEdges, 1) - 1, 100));
curDigits = ceil(log10(1 + numel(curIds)));

for curIdx = 1:numel(curIds)
    curId = curIds(curIdx);
    curEdge = mergedEdges(curId, 1:2);
    curScore = mergedEdges(curId, 3);
    
    curName = sprintf( ...
        '%0*d. %d-%d. Score: %.2f', ...
        curDigits, curIdx, curEdge, curScore);
    skel = skel.addTree(curName, box(:, 1)' + segPoints(curEdge, :) - 1);
end
