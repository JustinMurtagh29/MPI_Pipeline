function calculateGlobalCorrespondences(pSeg, cubeCoords1, cubeCoords2, segFile1, segFile2, bboxSmall1, bboxSmall2, saveFolder)
    % calculateGlobalCorrespondences -> Find correspondence between global segmentations by comparing over faces of adjacent cubes

    % Determine direction orthogonal to plane that touches in bboxSmall
    dirOfOverlap = find(~all(bboxSmall1 == bboxSmall2,2));

    % Load both segmentations (only face that touches)
    % This assumes cube1 is one with lower coordinates, see correspondenceFinderGlobal
    thisBbox1 = bboxSmall1;
    thisBbox1(dirOfOverlap,1) = thisBbox1(dirOfOverlap,2);
    touchFace1 = loadSegDataGlobal(pSeg, thisBbox1);
    thisBbox2 = bboxSmall2;
    thisBbox2(dirOfOverlap,2) = thisBbox2(dirOfOverlap,1);
    touchFace2 = loadSegDataGlobal(pSeg, thisBbox2);

    % create one dimensional array to compare both overlaps
    correspondences = cat(2, touchFace1(:), touchFace2(:));
    % remove correspondences to 0 (background)
    correspondences(any(correspondences == 0,2),:) = [];
    % do some statistics on correspondences
    [uniqueCorrespondences, ~, idxC] = unique(correspondences, 'rows');
    countsC = histc(idxC, 1:size(uniqueCorrespondences,1)); 
    [uniqueSegments, ~, idxS] = unique(correspondences(:)); 
    countsS = histc(idxS, 1:size(uniqueSegments,1)); 
    countsCnorm = countsC ./ max(arrayfun(@(x)countsS(uniqueSegments == x), uniqueCorrespondences),[],2);

    Util.save([saveFolder num2str(cubeCoords1, '%.2i') num2str(cubeCoords2, '%.2i') '.mat'], ...
        uniqueCorrespondences, countsC, countsCnorm, uniqueSegments, countsS);

end

