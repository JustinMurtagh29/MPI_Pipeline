function agglomerateDirectionality(axonsFinal, graph, segmentMeta, borderMeta, globalSegmentPCA)
    minSize = 2500;
    axonsFinalAll = [axonsFinal
        num2cell(setdiff(find(segmentMeta.axonProb > 0.5 & segmentMeta.voxelCount >= minSize), cell2mat(axonsFinal)))];
    disp('measuring agglo size');
    axonsFinalAllSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axonsFinalAll);
    axonsFinalAll(axonsFinalAllSize < minSize) = [];
    disp('preallocating direction score');
    directionScore = nan(length(borderMeta.borderSize), 2);
    selection = randperm(length(axonsFinalAll), 10);
    for idx = selection %1 : length(axonsFinalAll)
        idx
        currentAgglo = axonsFinalAll{idx};
        for idx2 = 1 : length(currentAgglo)
            % detect local surround
            surround = intersect(graph.neighbours{currentAgglo(idx2)}', currentAgglo);
            thisbbox = segmentMeta.box(:, :, currentAgglo(idx2));
            otherbboxes = segmentMeta.box(:, : , currentAgglo);

            surround = unique([currentAgglo(idx2); surround; currentAgglo(cellfun(@(x)bboxOverlap(thisbbox, x), num2cell(otherbboxes, [1,2])))]);
            % calculate PCA of local surround (Alessandro)
            massesIn = segmentMeta.voxelCount(surround);
            comVecsIn = bsxfun(@times, segmentMeta.centroid(:, surround)', [11.24, 11.24, 28]);
            covMatsIn = reshape(globalSegmentPCA.covMat(surround, :), [length(surround), 3, 3]);
            agglos = {1: length(surround)};
            [massesOut, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments(massesIn, comVecsIn, covMatsIn, agglos);
            [mypca, latent] = eig(covMatsOut);
            [latent1, idxLatent1] = max(sum(latent));
            if latent1 < 0.7
                continue;
            end
            % find all outgoing edges of current segment
            % outgoing = setdiff(graph.neighbours{currentAgglo(idx2)}(~isnan(graph.neighBordIdx{currentAgglo(idx2)})), currentAgglo);
            edgeIdxs = find(any(ismember(graph.edges,currentAgglo),2));
            borderIdxs = graph.borderIdx(edgeIdxs);
            outgoing = ~all(ismember(graph.edges(edgeIdxs, :), currentAgglo), 2) & ~isnan(borderIdxs);
            currentOutgoing = outgoing & any(ismember(graph.edges(edgeIdxs, :), currentAgglo(idx2)), 2);

            % calculate minmax score of those
            borderCoMs = bsxfun(@times, single(borderMeta.borderCoM(borderIdxs(outgoing), :)), [11.24, 11.24, 28]);

            borderCoMsLocalized = bsxfun(@minus, borderCoMs, comVecsOut);
            result = borderCoMsLocalized * mypca;
            scorePre = (result(:, idxLatent1) - min(result(:, idxLatent1))) / (max(result(:, idxLatent1)) - min(result(:, idxLatent1))) * 2 - 1;
            score = scorePre(currentOutgoing(outgoing));
            filename = ['/gaba/scratch/kboerg/direction/' num2str(idx,'%.5i') '_' num2str(idx2, '%.5i') '.nml'];
            nodesHere = {double(borderMeta.borderCoM(borderIdxs(currentOutgoing), :))};
            treenames =  {['tree1' num2str(idx,'%.5i') '_' num2str(idx2, '%.5i')]};
            comments = {arrayfun(@num2str,score,'uni', 0)};
            connectEM.generateSkeletonFromNodes(filename, nodesHere, treenames, comments);
            % directionScore(sub2ind(size(directionScore),graph.neighBordIdx(ismember(graph.neighbours, outgoing), (outgoing < currentAgglo(idx2))+1))) = score;
        end
    end
end

function overlap = bboxOverlap(bbox1, bbox2)
    threshold = 500;
    bbox1 = bsxfun(@times, bbox1', [11.24, 11.24, 28]);
    bbox2 = bsxfun(@times, bbox2', [11.24, 11.24, 28]);
    bbox1 = [bbox1(1, :) - threshold; bbox1(2, :) + threshold];
    bbox2 = [bbox2(1, :) - threshold; bbox2(2, :) + threshold];
    [A,B,C] = ndgrid(1:2,1:2,1:2);
    cornerIdx = cat(2, A(:), B(:), C(:)) +  repmat(0:2:4, 8, 1);
    corners = bbox1(cornerIdx);
    overlap = any(all(bsxfun(@gt, corners, bbox2(1,:)), 2) & all(bsxfun(@lt, corners, bbox2(2,:)), 2));
end
