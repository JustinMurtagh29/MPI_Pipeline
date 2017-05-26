function y = agglomerateDirectionality2(agglos, graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, voxelSize)

    % Preallocation
    y.latent = cell(numel(agglos),1); 
    y.pca = cell(numel(agglos),1);
    y.neighbours = cell(numel(agglos),1);
    y.prob = cell(numel(agglos),1);
    y.borderIdx = cell(numel(agglos),1);
    y.scores = cell(numel(agglos),1);

    % Loop over all agglomerates
    for idx1 = 1:length(agglos)

        % Remove segments from currentAgglo if covMatsIn contains NaN or Inf (small segments) 
        currentAgglo = sort(agglos{idx1});
        covMatsIn = globalSegmentPCA.covMat(currentAgglo, :);
        idx = any(isnan(covMatsIn) | isinf(covMatsIn),2);
        currentAgglo(idx) = [];
        clear covMatsIn idx;

        % Preallocate as empty as we do not know number of outgoing edges
        latent = [];
        pca = []; 
        neighbours = [];
        prob = [];
        borderIdx = [];
        scores = [];

        % If more than 2 segments in agglo -> find the surround, otherwise whole agglo is surround
        if numel(currentAgglo) > 2
            allbboxes = segmentMeta.box(:, : , currentAgglo);
            overlap = bboxOverlap(allbboxes, bboxDist, voxelSize');
            if all(overlap(:))
                surround = {currentAgglo};
            else
                directNeighbours = cellfun(@(x)ismember(currentAgglo, x), graph.neighbours(currentAgglo), 'uni', 0);
                directNeighbours = cat(2, directNeighbours{:});
                neighbourhood = (directNeighbours | overlap)';
                neighbourhood = unique(neighbourhood, 'rows');
                surround = cellfun(@(x)currentAgglo(x), mat2cell(neighbourhood, ones(size(neighbourhood,1),1), size(neighbourhood,2)), 'uni', 0);
                clear directNeighbours neighbourhood;
            end
            clear allbboxes overlap;
        else
            surround = {currentAgglo};
        end

        % Keep only nonempty surrounds
        idx = cellfun(@isempty, surround);
        surround(idx) = [];
        clear idx;

        % Combine PCAs for all surrounds in current agglo
        if numel(currentAgglo) > 1
            massesIn = segmentMeta.voxelCount(currentAgglo);
            comVecsIn = bsxfun(@times, segmentMeta.centroid(:, currentAgglo)', voxelSize);
            covMatsIn = reshape(globalSegmentPCA.covMat(currentAgglo, :), [length(currentAgglo), 3, 3]);
            agglos = cellfun(@(x)find(ismember(currentAgglo,x)), surround, 'uni', 0);
            [~, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments(massesIn, comVecsIn, covMatsIn, agglos);
            clear massesIn comVecsIn covMatsIn agglos;
        else
            comVecsOut = bsxfun(@times, segmentMeta.centroid(:, currentAgglo)', voxelSize);
            covMatsOut = reshape(globalSegmentPCA.covMat(currentAgglo, :), [1, 3, 3]);
        end

        % Calculate measures for each surround of current agglo
        for idx2=1:length(surround)

            % Calculate PCA from covariance matrix
            [thisPca, thisLatent] = pcaFromCovMat(squeeze(covMatsOut(idx2,:,:)));

            % Find all outgoing edges of current surround
            % without correspondences as these have no borderIdx and thereby CoM currently
            borderIdxs = cat(1, graph.neighBorderIdx{surround{idx2}});
            borderSegId = cat(1, graph.neighbours{surround{idx2}});
            borderProb = cat(1, graph.neighProb{surround{idx2}});
            % Indices to border* that are not correspondences
            outgoing = ~isnan(borderIdxs) & ~ismember(borderSegId, surround{idx2});
            % Indices to border* that originate for current segment
            currentOutgoing = outgoing & ~ismember(borderSegId, currentAgglo);

            if ~any(currentOutgoing)
                continue;
            end

            % Get center of mass of all border of surround
            borderCoMs = bsxfun(@times, single(borderMeta.borderCoM(borderIdxs(outgoing), :)), voxelSize);
            % Center on CoM of current agglo
            borderCoMsLocalized = bsxfun(@minus, borderCoMs, comVecsOut(idx2,:));
            % Project into PCA space of current agglo
            result = borderCoMsLocalized * thisPca;
            scorePre = (result(:,1) - min(result(:,1))) / (max(result(:,1)) - min(result(:,1))) * 2 - 1;
            score = scorePre(currentOutgoing(outgoing));

            % Collect output
            latent = [latent; repmat(thisLatent', size(score,1), 1)];
            pca = cat(3, pca, repmat(thisPca, 1, 1, size(score,1)));
            neighbours = [neighbours; borderSegId(currentOutgoing)];
            prob = [prob; borderProb(currentOutgoing)];
            borderIdx = [borderIdx; borderIdxs(currentOutgoing)];
            scores = [scores; score];

        end

        % Use (possible) redundancy in surround calculations to find "global" endings
        if numel(surround) > 1
            % We like to extract the minimal absolute score and in case of quality: highest latent score
            [~, sortIdx] = sortrows([abs(scores) latent(:,1)], [1 -2]);
            [~, idx] = unique(borderIdx(sortIdx));
            latent = latent(sortIdx(idx), :);
            pca = pca(:,:,sortIdx(idx));
            neighbours = neighbours(sortIdx(idx));
            prob = prob(sortIdx(idx));
            borderIdx = borderIdx(sortIdx(idx));
            scores = scores(sortIdx(idx));
            clear sortIdx idx;
        end

        % More collecting
        y.latent{idx1} = latent;
        y.pca{idx1} = pca;
        y.neighbours{idx1} = neighbours;
        y.prob{idx1} = prob;
        y.borderIdx{idx1} = borderIdx;
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
