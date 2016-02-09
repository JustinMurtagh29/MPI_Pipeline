function calculateLocalCorrespondences(cubeCoords1, cubeCoords2, segFile1, segFile2, bboxSmall1, bboxSmall2, bboxBig1, bboxBig2, nrPlanes, saveFolder )
    % calculateLocalCorrespondences -> Find correspondence between local segmentations by comparing overlap around face between cubes

    % Determine direction orthogonal to plane that touches in bboxSmall
    dirOfOverlap = find(~all(bboxSmall1 == bboxSmall2,2));

    % This assumes that bbox2 is larger coordinate one, see getOverlaps, also assumes non overlapping small bounding boxes
    planesNeeded = (bboxSmall1(dirOfOverlap,2)-nrPlanes+1):1:(bboxSmall2(dirOfOverlap,1)-nrPlanes+1);

    % Load both segmentations (bboxSmall region only + region needed for
    % overlap calculations)
    overlap1 = extractBboxSmall(segFile1, bboxSmall1, bboxBig1, dirOfOverlap, planesNeeded);
    overlap2 = extractBboxSmall(segFile2, bboxSmall2, bboxBig2, dirOfOverlap, planesNeeded);

    % Make dimension orthogonal to plane last dimension in array
    permuteVector = [1 2 3];
    permuteVector(permuteVector == dirOfOverlap) = 3;
    permuteVector(3) = dirOfOverlap;
    overlap1 = permute(overlap1, permuteVector);
    overlap2 = permute(overlap2, permuteVector);

    % Determine which ids are present in respecitve cube
    idsInSmall1 = unique(overlap1(:,:,1:nrPlanes));
    idsInSmall2 = unique(overlap2(:,:,end-nrPlanes+1:end));

    % create one dimensional array to compare both overlaps
    correspondences = [reshape(overlap1, [numel(overlap1) 1 1]) reshape(overlap2, [numel(overlap2) 1 1])];
    % set all ids to zeros that are not within bboxSmall of respective cube
    idxToSetZero = ~ismember(correspondences(:,1),idsInSmall1);
    correspondences(idxToSetZero,1) = 0;
    idxToSetZero = ~ismember(correspondences(:,2),idsInSmall2);
    correspondences(idxToSetZero,2) = 0;

    % remove correspondences to 0 (background or set above)
    correspondences(any(correspondences == 0,2),:) = [];

    % find unique correspondences
    uniqueCorrespondences = unique(correspondences, 'rows');

    % Create main correspondence folder if necessary
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    save([saveFolder num2str(cubeCoords1, '%.2i') num2str(cubeCoords2, '%.2i') '.mat'], 'cubeCoords1', 'cubeCoords2', 'uniqueCorrespondences');

end

