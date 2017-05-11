function y = agglomerateDirectionality(axonsFinalAll, graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, visualize)
    y.latent = sparse(1, max(cell2mat(axonsFinalAll)));
    y.edges = [];
    y.scores = [];
    axonsFinalAllSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axonsFinalAll);
    minSize = 100;
    axonsFinalAll(axonsFinalAllSize < minSize) = [];
    for idx = 1 : length(axonsFinalAll)
        currentAgglo = axonsFinalAll{idx};
        for idx2 = 1 : length(currentAgglo)
            % detect local surround
            surround = intersect(graph.neighbours{currentAgglo(idx2)}', currentAgglo);
            thisbbox = segmentMeta.box(:, :, currentAgglo(idx2));
            otherbboxes = segmentMeta.box(:, : , currentAgglo);

            surround = unique([currentAgglo(idx2); surround; currentAgglo(cellfun(@(x)bboxOverlap(thisbbox, x, bboxDist), num2cell(otherbboxes, [1,2])))]);
            covMatsIn = num2cell(reshape(globalSegmentPCA.covMat(surround, :), [length(surround), 3, 3]), [2 3]);
            surround = surround(~cellfun(@(x)any(isnan(x(:)) | isinf(x(:))), covMatsIn));
            if isempty(surround)
                continue
            end
            % calculate PCA of local surround (Alessandro)
            massesIn = segmentMeta.voxelCount(surround);
            comVecsIn = bsxfun(@times, segmentMeta.centroid(:, surround)', [11.24, 11.24, 28]);
            covMatsIn = reshape(globalSegmentPCA.covMat(surround, :), [length(surround), 3, 3]);
            agglos = {1: length(surround)};
            [massesOut, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments(massesIn, comVecsIn, covMatsIn, agglos);
            [mypca, latent] = eig(squeeze(covMatsOut));
            assert(min(latent(:))>-1E5);
            latent = latent / sum(sum(latent));
            [latent1, idxLatent1] = max(sum(latent));

            % find all outgoing edges of current surround
            borderIdxs = cat(1, graph.neighBorderIdx{surround});
            borderSegId = cat(2, graph.neighbours{surround});
            borderLookUp = repelem(surround', cellfun(@length, graph.neighbours(surround)))';
            outgoing = ~isnan(borderIdxs) & ~ismember(borderSegId', surround);
            currentOutgoing = outgoing & borderLookUp == currentAgglo(idx2);
            if ~any(currentOutgoing)
                continue;
            end

            % calculate minmax score of those
            borderCoMs = bsxfun(@times, single(borderMeta.borderCoM(borderIdxs(outgoing), :)), [11.24, 11.24, 28]);

            borderCoMsLocalized = bsxfun(@minus, borderCoMs, comVecsOut);
            result = borderCoMsLocalized * mypca;
            scorePre = (result(:, idxLatent1) - min(result(:, idxLatent1))) / (max(result(:, idxLatent1)) - min(result(:, idxLatent1))) * 2 - 1;
            score = scorePre(currentOutgoing(outgoing));
            y.latent(currentAgglo(idx2)) = latent1;
            y.edges = [y.edges; repmat(currentAgglo(idx2), size(score)), borderSegId(currentOutgoing)'];
            y.scores = [y.scores; score];
            if visualize
                borderProb = borderIdxs; %cat(1, graph.neighProb{surround});
                currentBorderProb = borderProb(currentOutgoing);
                currentBorderIdxs = borderIdxs(currentOutgoing);
                treename= ['size' num2str(sum(segmentMeta.voxelCount(surround))) '_latent' num2str(latent1)];
                if latent1 < 0.7
                    treename = [treename, 'unused'];
                end

                filename = ['/gaba/scratch/kboerg/direction/' num2str(idx,'%.5i') '_' num2str(idx2, '%.5i') '.nml'];
                nodesHere = {double(borderMeta.borderCoM(borderIdxs(currentOutgoing), :))};
                treenames =  {[treename 'tree1' num2str(idx,'%.5i') '_' num2str(idx2, '%.5i')]};

                comments = {arrayfun(@(x)['score_' num2str(score(x)), '_p_' num2str(currentBorderProb(x)) '_size_' num2str(borderMeta.borderSize(currentBorderIdxs(x)))], 1:length(score),'uni', 0)};
                connectEM.generateSkeletonFromNodes(filename, nodesHere, treenames, comments);
            end
        end
    end
end

function overlap = bboxOverlap(bbox1, bbox2, bboxDist)
    persistent cornerIdx;
    bbox1 = bsxfun(@times, bbox1', [11.24, 11.24, 28]);
    bbox2 = bsxfun(@times, bbox2', [11.24, 11.24, 28]);
    bbox1 = [bbox1(1, :) - bboxDist / 2; bbox1(2, :) + bboxDist / 2];
    bbox2 = [bbox2(1, :) - bboxDist / 2; bbox2(2, :) + bboxDist / 2];
    if isempty(cornerIdx)
        [A,B,C] = ndgrid(1:2,1:2,1:2);
        cornerIdx = cat(2, A(:), B(:), C(:)) +  repmat(0:2:4, 8, 1);
    end
    corners = bbox1(cornerIdx);
    overlap = any(all(bsxfun(@gt, corners, bbox2(1,:)), 2) & all(bsxfun(@lt, corners, bbox2(2,:)), 2));
end
