function y = agglomerateDirectionality(axonsFinalAll, graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, visualize)

    % Preallocation
    y.latent = cell(numel(axonsFinalAll),1);
    y.edges = cell(numel(axonsFinalAll),1);
    y.scores = cell(numel(axonsFinalAll),1);

    % Loop over all agglomerates and all segments within each agglomerate
    for idx1 = 1:length(axonsFinalAll)
        currentAgglo = axonsFinalAll{idx1};
        edges = [];
        scores = [];
        for idx2 = 1:length(currentAgglo)
            
            % If more than 2 segments in agglo -> find the surround, otherwise whole agglo is surround
            if numel(currentAgglo) > 2
                % Detect local surround: All direct neighbours
                surround = intersect(graph.neighbours{currentAgglo(idx2)}', currentAgglo);
                thisbbox = segmentMeta.box(:, :, currentAgglo(idx2));
                otherbboxes = segmentMeta.box(:, : , currentAgglo);
                % Add to local surround: All segments whose bouding box overlap when each bounding box is expanded
                closeSegments = currentAgglo(bboxOverlap(thisbbox, otherbboxes, bboxDist));
                surround = unique([currentAgglo(idx2); surround; closeSegments]);
            else
                surround = currentAgglo;    
            end

            % Remove segments from surround if covMatsIn is NaN (small segments)
            covMatsIn = globalSegmentPCA.covMat(surround, :);
            surround = surround(any(isnan(covMatsIn) | isinf(covMatsIn),2));
            % Skip current segment if surround is empty
            if isempty(surround)
                continue
            end

            % Calculate PCA of local surround
            if numel(surround) > 1
                massesIn = segmentMeta.voxelCount(surround);
                comVecsIn = bsxfun(@times, segmentMeta.centroid(:, surround)', [11.24, 11.24, 28]);
                covMatsIn = reshape(globalSegmentPCA.covMat(surround, :), [length(surround), 3, 3]);
                agglos = {1: length(surround)};
                [~, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments(massesIn, comVecsIn, covMatsIn, agglos);
            else
                comVecsOut = bsxfun(@times, segmentMeta.centroid(:, surround)', [11.24, 11.24, 28]);
                covMatsOut = squeeze(reshape(globalSegmentPCA.covMat(surround, :), [length(surround), 3, 3]));
            end
            [mypca, latent] = eig(covMatsOut);
            assert(min(latent(:))>-1E5);
            latent = latent / sum(sum(latent));
            [latent1, idxLatent1] = max(sum(latent));

            % Find all outgoing edges of current segment
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

            % Save output in output structure
            y.latent{idx1}(idx2) = latent1;
            edges = [edges; repmat(currentAgglo(idx2), size(score)), borderSegId(currentOutgoing)'];
            scores = [scores; score];

            % Write out skeletons if flag is set
            if visualize
                borderProb = cat(1, graph.neighProb{surround});
                currentBorderProb = borderProb(currentOutgoing);
                currentBorderIdxs = borderIdxs(currentOutgoing);
                treename= ['size' num2str(sum(segmentMeta.voxelCount(surround))) '_latent' num2str(latent1)];
                if latent1 < 0.7
                    treename = [treename, 'unused'];
                end
                filename = ['/gaba/scratch/kboerg/direction/' num2str(idx,'%.5i') '_' num2str(idx2, '%.5i') '.nml'];
                nodesHere = {double(borderMeta.borderCoM(borderIdxs(currentOutgoing), :))};
                treenames =  {[treename 'tree1' num2str(idx,'%.5i') '_' num2str(idx2, '%.5i')]};
                comments = {arrayfun(@(x)['score_' num2str(score(x)), '_p_' num2str(currentBorderProb(x)) ...
                    '_size_' num2str(borderMeta.borderSize(currentBorderIdxs(x)))], 1:length(score),'uni', 0)};
                connectEM.generateSkeletonFromNodes(filename, nodesHere, treenames, comments);
            end
        end
        y.edges{idx1} = edges;
        y.scores{idx1} = scores;
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

