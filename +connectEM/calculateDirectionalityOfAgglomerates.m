function y = calculateDirectionalityOfAgglomerates(agglos, graph, segmentMeta, borderMeta, globalSegmentPCA, options)
% Calculate directionality of agglomerates
% Basic idea is to always look at a local surround of each segment in the agglomerate 
% and calculate its principal component and a score that shows for each border
% pointing out of this agglomerate whether it is close to an ending

    % Preallocation for storing results for each agglomerate
    y.neighbours = cell(numel(agglos),1);
    y.prob = cell(numel(agglos),1);
    y.borderIdx = cell(numel(agglos),1);
    y.latent = cell(numel(agglos),1); 
    y.pca = cell(numel(agglos),1);
    y.scores = cell(numel(agglos),1);

    % Loop over all agglomerates
    for idx1 = 1:length(agglos)

        % Remove segments from currentAgglo if covMatsIn contains NaN or Inf (small segments) 
        currentAgglo = sort(agglos{idx1});
        covMatsIn = globalSegmentPCA.covMat(currentAgglo, :);
        idx = any(isnan(covMatsIn) | isinf(covMatsIn),2);
        currentAgglo(idx) = [];
        covMatsIn(idx,:) = [];
        covMatsIn = reshape(covMatsIn, [length(currentAgglo), 3, 3]);
        clear idx;

        % If more than 2 segments in agglo -> find the surround, otherwise whole agglo is surround
        if numel(currentAgglo) > 2
            allbboxes = segmentMeta.box(:, : , currentAgglo);
            overlap = bboxOverlap(allbboxes, options.bboxDist, options.voxelSize');
            if all(overlap(:))
                surround = {currentAgglo};
                neighbourhood = {true(1,numel(currentAgglo))};
            else
                directNeighbours = cellfun(@(x)ismember(currentAgglo, x), graph.neighbours(currentAgglo), 'uni', 0);
                directNeighbours = cat(2, directNeighbours{:});
                neighbourhood = (directNeighbours | overlap)';
                neighbourhood = unique(neighbourhood, 'rows');
                neighbourhood = mat2cell(neighbourhood, ones(size(neighbourhood,1),1), size(neighbourhood,2));
                surround = cellfun(@(x)currentAgglo(x), neighbourhood, 'uni', 0);
                clear directNeighbours;
            end
            clear allbboxes overlap;
        else
            surround = {currentAgglo};
            neighbourhood = {true(1,numel(currentAgglo))};
        end

        % Combine PCAs for all surrounds in current agglo
        if numel(currentAgglo) > 1
            massesIn = segmentMeta.voxelCount(currentAgglo);
            comVecsIn = bsxfun(@times, segmentMeta.centroid(:, currentAgglo)', options.voxelSize);
            surroundLocal = cellfun(@(x)find(x), neighbourhood, 'uni', 0);
            [~, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments(massesIn, comVecsIn, covMatsIn, surroundLocal);
            clear massesIn comVecsIn;
        else
            comVecsOut = bsxfun(@times, segmentMeta.centroid(:, currentAgglo)', options.voxelSize);
            covMatsOut = reshape(globalSegmentPCA.covMat(currentAgglo, :), [1, 3, 3]);
        end
        clear covMatsIn;
        
        % Get representation of all borders out of current agglomerate
        borderSegId = cat(1, graph.neighbours{currentAgglo});
        borderProb = cat(1, graph.neighProb{currentAgglo}); 
        borderIdxs = cat(1, graph.neighBorderIdx{currentAgglo});
        % Limit to borders that are not correspondences and point outside current agglomerate
        outgoing = ~isnan(borderIdxs) & ~ismember(borderSegId, currentAgglo);
        borderSegId = borderSegId(outgoing);
        borderProb = borderProb(outgoing); 
        borderIdxs = borderIdxs(outgoing);
        % Number of outgoing non correspondence borders
        nrBorder = numel(borderIdxs);
        clear outgoing;

        % Preallocate arrays for this agglomerate
        latent = zeros(nrBorder,3);
        pca = zeros(3,3,nrBorder); 
        scores = zeros(nrBorder,1);

        % Calculate measures for each surround of current agglo
        for idx2=1:length(surround)

            % Calculate PCA from covariance matrix
            [thisPca, thisLatent] = pcaFromCovMat(squeeze(covMatsOut(idx2,:,:)));

            % Get center of mass of all border of surround
            borderCoMs = bsxfun(@times, single(borderMeta.borderCoM(borderIdxs(surroundLocal{idx2}), :)), options.voxelSize);
            % Center on CoM of current surround
            borderCoMsLocalized = bsxfun(@minus, borderCoMs, comVecsOut(idx2,:));
            % Project into PCA space of current surround
            result = borderCoMsLocalized * thisPca;
            score = (result(:,1) - min(result(:,1))) / (max(result(:,1)) - min(result(:,1))) * 2 - 1;

            % Collect output
            idx = abs(scores(surroundLocal{idx2})) < abs(score);
            nrReplace = sum(idx);
            latent(surroundLocal{idx2}(idx),:) = repmat(thisLatent', nrReplace, 1);
            pca(:,:,surroundLocal{idx2}(idx)) = repmat(thisPca, 1, 1, nrReplace);
            scores(surroundLocal{idx2}(idx)) = score(idx);

        end

        % More collecting
        y.neighbours{idx1} = borderSegId;
        y.prob{idx1} = borderProb;
        y.borderIdx{idx1} = borderIdxs;
        y.latent{idx1} = latent;
        y.pca{idx1} = pca;
        y.scores{idx1} = scores;

    end
end

function overlap = bboxOverlap(bbox1, bboxDist, voxelSize)
% Get symetric overlap matrix for a given set of bounding boxes
% Will include corner, edges and face only overlaps

    bbox1 = bsxfun(@times, bbox1, voxelSize);
    bbox1 = bsxfun(@plus, bbox1, [-bboxDist/2, bboxDist/2]);
    bbox2 = permute(bbox1, [1 2 4 3]);
    forwDir = any(bsxfun(@gt, bbox1(:,1,:), bbox2(:,2,:,:)),1);
    backDir = any(bsxfun(@lt, bbox1(:,2,:), bbox2(:,1,:,:)),1);
    overlap = ~squeeze(forwDir | backDir);

end

function [pca, latent] = pcaFromCovMat(covMat)
% Calculate PCA and cumulative energy from covariance matrix

    [pca, latent] = eig(covMat);
    latent = diag(latent);
    latent = latent ./ sum(latent);
    [latent, idx] = sort(latent, 'descend');
    pca = pca(:,idx);

end

