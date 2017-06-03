function q = generateAxonQueries(p, graph, segmentMeta, borderMeta, directionality, axons, options)

    if ~exist('options', 'var')
        options.latentScore = 0.7;
        options.segDirScore = 0.9;
        options.border = [3000; -3000];
        options.writeTasksToFile = true;
        options.boundingBoxForTasks = false;
    end

    % Find all endings (given criteria on latent, directionality, orientation along z-axis)
    correctlyFlippedScores = cellfun(@(x,y)x.*squeeze(sign(y(3,1,:))), directionality.scores, directionality.pca, 'uni', 0);
    idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, directionality.latent, 'uni', 0);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, directionality.scores, 'uni', 0);
    idxCorrectOrientation = cellfun(@(x)x > 0, correctlyFlippedScores, 'uni', 0);
    idxAll = cellfun(@(x,y,z)find(x&y&z), idxDirectional, idxEnding, idxCorrectOrientation, 'uni', 0);
    % Keep only largest score in each agglomerate for now
    nrCanidates = cellfun(@numel, idxAll);
    idxCanidateFound = nrCanidates > 0;
    [~, sortIdx] = cellfun(@(x,y)sort(x(y), 'descend'), correctlyFlippedScores(idxCanidateFound), idxAll(idxCanidateFound), 'uni', 0);
    idxUse = cellfun(@(x,y)x(y(1)), idxAll(idxCanidateFound), sortIdx, 'uni', 0);

    % Extract some values for borders we want to query
    thisPca = cell2mat(cellfun(@(x,y)x(:,1,y)', directionality.pca(idxCanidateFound), idxUse, 'uni', 0));
    thisScores = cellfun(@(x,y)x(y), directionality.scores(idxCanidateFound), idxUse);
    thisBorderIdx = cellfun(@(x,y)x(y), directionality.borderIdx(idxCanidateFound), idxUse);
    axons = axons(idxCanidateFound);

    % check for border of dataset
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm));
    bboxSmall = p.bbox + borderVoxel';
    borderPositions = double(borderMeta.borderCoM(thisBorderIdx,:));
    outsideBbox = ~(all(bsxfun(@gt, borderPositions, bboxSmall(:, 1)'), 2) & ...
        all(bsxfun(@lt, borderPositions, bboxSmall(:, 2)'), 2));

    % Position and direction of queries
    borderPositions = mat2cell(borderPositions(~outsideBbox,:), ones(sum(~outsideBbox),1), 3);
    directions = bsxfun(@times, bsxfun(@times, sign(thisScores(~outsideBbox)), thisPca(~outsideBbox,:)), 1 ./ p.raw.voxelSize);
    directions = mat2cell(directions, ones(size(directions,1),1), 3);
    axons = axons(~outsideBbox);
    tic
    positions = cellfun(@(x,y,z)connectEM.correctQueryLocationToEndOfSegment(p, x, y, z, 200), axons, borderPositions, directions, 'uni', 0);
    toc

    % Calculate euler angles & put into old format
    q.pos = positions;
    q.dir = directions;
    tic
    [phi, theta, psi] = cellfun(@(x)calculateEulerAngles(x, p.raw.voxelSize), q.dir);
    toc
    q.angles = mat2cell(cat(2, phi, theta, psi), ones(numel(phi),1), 3);

end

