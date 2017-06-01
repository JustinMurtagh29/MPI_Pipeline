function generateAxonQueries(p, graph, segmentMeta, borderMeta, directionality, axons, options)

    if ~exist('options', 'var')
        options.latentScore = 0.7;
        options.segDirScore = 0.9;
        options.border = [-3000 3000];
        options.writeTasksToFile = true;
        options.boundingBoxForTasks = false;
    end

    % Find all endings (given criteria on latent, directionality, orientation along z-axis)
    correctlyFlippedScores = cellfun(@(x,y)x.*squeeze(sign(y(3,1,:))), directionality.scores, directionality.pca, 'uni', 0);
    idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, directionality.latent, 'uni', 0);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, directionality.scores, 'uni', 0);
    idxCorrectOrientation = cellfun(@(x)x > 0, correctlyFlippedScores, 'uni', 0);
    idxAll = cellfun(@(x,y,z)find(x&y&z), idxDirectional, idxEnding, idxCorrectOrientation, 'uni', 0);
    % Keep only largest score in each agglomerate
    [~, sortIdx] = cellfun(@(x,y)sort((x(y), 'descend')), correctlyFlippedScores, idxAll, 'uni', 0);
    idxUse = cellfun(@(x,y)x(y(1)), idxAll, sortIdx, 'uni', 0);

    % Extract some values for borders we want to query
    thisPca = cellfun(@(x,y)x(:,1,y), directionality.pca, idxUse);
    thisScores = cellfun(@(x,y)x(y), directionality.scores, idxUse);
    thisBorderIdx = cellfun(@(x,y)x(y), directionality.borderIdx, idxUse);

    % Position and direction of queries
    directions = sign(thisScores) .* thisPca ./ p.raw.voxelSize;
    positions = cellfun(@(x,y)connectEM.correctQueryLocationToEndOfSegment(p, x, y, direction, 200), ...
        axons, double(borderMeta.borderCoM(borderIdx)), direction, 200, struct('failGracefully', true));

    % check for border of dataset
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm);
    bbox = p.bbox + borderVoxel;
    outsideBbox = ~(all(bsxfun(@gt, direction, bboxSmall(:, 1)'), 2) & ...
                    all(bsxfun(@lt, direction, bboxSmall(:, 2)'), 2));
    positions(outsideBBox,:) = [];
    directions(outsideBBox,:) = [];

    % Write queries to file
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    connectEM.generateQueriesFromData(p, segmentMeta, q, outputFolder, options)

end

