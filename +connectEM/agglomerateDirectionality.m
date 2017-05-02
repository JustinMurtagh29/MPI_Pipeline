function agglomerateDirectionality(axonsFinal, graph, segmentMeta, borderMeta)
axonsFinalAll = [axonsFinal ...
    num2cell(setdiff(find(segmentMeta.axonProb > 0.5), cell2mat(axonsFinal)))];
axonFinalAllSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axonsFinalAll);
axonsFinalAll(axonsFinalAllSize < 2500) = [];
directionScore = nan(length(borderMeta.borderSize), 2);
for idx = 1 : length(axonsFinalAll)
    currentAgglo = axonsFinalAll{idx};
    for idx2 = 1 : length(currentAgglo)
        % detect local surround
        surround = intersect(graph.neighbours{currentAgglo(idx2)}, currentAgglo);
        thisbbox = segmentMeta.box(:, :, currentAgglo(idx2));
        otherbboxes = segmentMeta.box(:, : , currentAgglo);

        surround = union(currentAgglo(idx2), surround, currentAgglo(bsxfun(@(x,y)overlap(x), otherbboxes, otherbboxes)));
        % calculate PCA of local surround (Alessandro)
        [pca, latent] = blabla;
        if latent(1) < 0.7
            continue;
        end
        % find all outgoing edges of current segment
        outgoing = setdiff(graph.neighbours{currentAgglo(idx2)}(~isnan(graph.neighBordIdx{currentAgglo(idx2)})), currentAgglo);
        % calculate minmax score of those
        localSurroundCoM = sum(bsxfun(@times, segmentMeta.point(:, surround), segmentMeta.voxelCount(surround))) / sum(segmentMeta.voxelCount(surround));
        borderCoMs = borderMeta.borderCoM(graph.neighBordIdx(ismember(graph.neighbours, outgoing)));
        borderCoMsLocalized = borderCoMs - localSurroundCoM;
        result = pca * borderComLocalized;
        score = (result(:, 1) - min(result(:, 1))) / (max(result(:, 1)) - min(result(:, 1))) * 2 - 1;

        directionScore(sub2ind(size(directionScore),graph.neighBordIdx(ismember(graph.neighbours, outgoing), (outgoing < currentAgglo(idx2))+1))) = score;
    end
end
end
function overlap = bboxOverlap(bbox1, bbox2)
    threshold = 500;
    bbox1 = bsxfun(@times, bbox1', [11.24, 11.24, 28]);
    bbox2 = bsxfun(@times, bbox2', [11.24, 11.24, 28]);
    bbox1 = [bbox1(1, :) - threshold; bbox1(2, :) + threshold];
    bbox2 = [bbox2(1, :) - threshold; bbox2(2, :) + threshold];
    [A,B,C] = ndgrid(1:2,1:2,1:2);
    cornerIdea =cat(2, A(:), B(:), C(:)) +  repmat([0:2:4], 8,1);
    corners = bbox1(cornerIdea);
    overlap = any(all(bsxfun(@gt, corners, bbox2(1,:)), 2) & all(bsxfun(@lt, corners, bbox2(2,:)), 2));
end
