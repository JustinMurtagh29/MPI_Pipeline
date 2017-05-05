function y = agglomerateDirectionality2(axonsFinalAll, graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist)

    % Preallocation
    y.latent = cell(numel(axonsFinalAll),1);
    y.edges = cell(numel(axonsFinalAll),1);
    y.scores = cell(numel(axonsFinalAll),1);
    y.segIds = cell(numel(axonsFinalAll),1);
    % Loop over all agglomerates and all segments within each agglomerate
    for idx1 = 1:length(axonsFinalAll)
        currentAgglo = axonsFinalAll{idx1};
        % Remove segments from currentAgglo if covMatsIn contains NaN or Inf (small segments)
        covMatsIn = globalSegmentPCA.covMat(currentAgglo, :);
        currentAgglo(any(isnan(covMatsIn) | isinf(covMatsIn),2)) = [];
        % Preallocation
        edges = [];
        scores = [];
        segIds = [];
        surround = cell(numel(currentAgglo), 1);
        for idx2 = 1:length(currentAgglo)
            % If more than 2 segments in agglo -> find the surround, otherwise whole agglo is surround
            if numel(currentAgglo) > 2
                % Detect local surround: All direct neighbours
                surround{idx2} = intersect(graph.neighbours{currentAgglo(idx2)}', currentAgglo);
                thisbbox = segmentMeta.box(:, :, currentAgglo(idx2));
                otherbboxes = segmentMeta.box(:, : , currentAgglo);
                % Add to local surround: All segments whose bouding box overlap when each bounding box is expanded
                closeSegments = currentAgglo(bboxOverlap(thisbbox, otherbboxes, bboxDist));
                surround{idx2} = unique([currentAgglo(idx2); surround{idx2}; closeSegments]);
            else
                surround{idx2} = currentAgglo;
            end
        end
        % Keep only nonempty surrounds
        idx = cellfun(@isempty, surround);
        surround(idx) = [];
        currentAgglo(idx) = [];
        if numel(currentAgglo > 1)
            massesIn = segmentMeta.voxelCount(currentAgglo);
            comVecsIn = bsxfun(@times, segmentMeta.centroid(:, currentAgglo)', [11.24, 11.24, 28]);
            covMatsIn = reshape(globalSegmentPCA.covMat(currentAgglo, :), [length(currentAgglo), 3, 3]);
            agglos = cellfun(@(x)find(ismember(currentAgglo,x)), surround, 'uni', 0);
            [~, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments(massesIn, comVecsIn, covMatsIn, agglos);
        else
            comVecsOut = bsxfun(@times, segmentMeta.centroid(:, surround)', [11.24, 11.24, 28]);
            covMatsOut = squeeze(reshape(globalSegmentPCA.covMat(surround, :), [length(surround), 3, 3]));
        end
        for idx2=1:length(surround)
            [mypca, latent] = eig(squeeze(covMatsOut(idx2,:,:)));
            assert(min(latent(:))>-1E5);
            latent = latent / sum(sum(latent));
            [latent1, idxLatent1] = max(sum(latent));

            % Find all outgoing edges of current segment
            borderIdxs = cat(1, graph.neighBorderIdx{surround{idx2}});
            borderSegId = cat(2, graph.neighbours{surround{idx2}});
            borderLookUp = repelem(surround{idx2}', cellfun(@length, graph.neighbours(surround{idx2})))';
            outgoing = ~isnan(borderIdxs) & ~ismember(borderSegId', surround{idx2});
            currentOutgoing = outgoing & borderLookUp == currentAgglo(idx2);
            if ~any(currentOutgoing)
                continue;
            end

            % calculate minmax score of those
            borderCoMs = bsxfun(@times, single(borderMeta.borderCoM(borderIdxs(outgoing), :)), [11.24, 11.24, 28]);
            borderCoMsLocalized = bsxfun(@minus, borderCoMs, comVecsOut(idx2,:));
            result = borderCoMsLocalized * mypca;
            scorePre = (result(:, idxLatent1) - min(result(:, idxLatent1))) / (max(result(:, idxLatent1)) - min(result(:, idxLatent1))) * 2 - 1;
            score = scorePre(currentOutgoing(outgoing));

            % Collect output
            y.latent{idx1}(idx2) = latent1;
            edges = [edges; repmat(currentAgglo(idx2), size(score)), borderSegId(currentOutgoing)'];
            scores = [scores; score];
            segIds = [segIds; currentAgglo(idx2)];
        end
        y.edges{idx1} = edges;
        y.scores{idx1} = scores;
        y.segIds{idx1} = segIds;
    end
end

function overlap = bboxOverlap(bbox1, bbox2, bboxDist)
    persistent cornerIdx;
    bbox1 = bsxfun(@times, bbox1, [11.24; 11.24; 28]);
    bbox2 = bsxfun(@times, bbox2, [11.24; 11.24; 28]);
    bbox1 = bsxfun(@plus, bbox1, [-bboxDist/2, bboxDist/2]);
    bbox2 = bsxfun(@plus, bbox2, [-bboxDist/2, bboxDist/2]);
    if isempty(cornerIdx)
        [A,B,C] = ndgrid([1 4],[2 5],[3 6]);
        cornerIdx = cat(2, A(:), B(:), C(:));
    end
    corners = bbox1(cornerIdx)';
    belowUpper = all(bsxfun(@lt, corners, bbox2(:,2,:)),1);
    aboveLower = all(bsxfun(@gt, corners, bbox2(:,1,:)),1);
    within = any(belowUpper & aboveLower,2);
    overlap = squeeze(within);
end

