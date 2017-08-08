function calculateGlobalCorrespondences(pSeg, cubeCoords1, cubeCoords2, bboxSmall1, bboxSmall2, saveFolder)
    % calculateGlobalCorrespondences -> Find correspondence between global segmentations by comparing over faces of adjacent cubes

    % Determine direction orthogonal to plane that touches in bboxSmall
    dirOfOverlap = find(~all(bboxSmall1 == bboxSmall2,2));

    % Load both segmentations (only face that touches)
    % This assumes cube1 is one with lower coordinates, see correspondenceFinderGlobal
    bboxSmall1(dirOfOverlap,1) = bboxSmall1(dirOfOverlap,2);
    touchFace1 = loadSegDataGlobal(pSeg, bboxSmall1);
    bboxSmall2(dirOfOverlap,2) = bboxSmall2(dirOfOverlap,1);
    touchFace2 = loadSegDataGlobal(pSeg, bboxSmall2);

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
    com = getCorrespondenceCom(touchFace1, touchFace2, ...
        uniqueCorrespondences, bboxSmall1);
    
    Util.save([saveFolder num2str(cubeCoords1, '%.2i') num2str(cubeCoords2, '%.2i') 'global.mat'], ...
        uniqueCorrespondences, countsC, countsCnorm, uniqueSegments, countsS, com);
end

function coms = getCorrespondenceCom(touchFace1, touchFace2, ...
    uniqueCorrespondences, bboxSmall1)
% Get correspondence com/points by taking the mean of the point with
% smallest bwdist of walls for both correspondence segments.
% Very rudimentary implementation that could possibly be speed up.
% OUTPUT coms: [Nx3] int
%           Global coms of the correpondence.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

[face1, ids1] = Seg.Local.fromGlobal(touchFace1);
global2Local1 = sparse(double(ids1), 1, 1:length(ids1));
[face2, ids2] = Seg.Local.fromGlobal(touchFace2);
global2Local2 = sparse(double(ids2), 1, 1:length(ids2));
stats1 = regionprops(face1, 'PixelIdxList');
stats2 = regionprops(face2, 'PixelIdxList');
bw1 = bwdist(touchFace1 == 0);
bw2 = bwdist(touchFace2 == 0);
coms = zeros(size(uniqueCorrespondences, 1), 3);
for i = 1:size(uniqueCorrespondences, 1)
    plist1 = stats1(global2Local1(uniqueCorrespondences(i, 1))).PixelIdxList;
    plist2 = stats2(global2Local2(uniqueCorrespondences(i, 2))).PixelIdxList;
    % operate on intersection of the pixels only
    plist = intersect(plist1, plist2); 
    [~,idx1] = max(bw1(plist));
    p1 = plist(idx1);
    [x1, y1, z1] = ind2sub(size(touchFace1), p1);
    [~,idx2] = max(bw2(plist));
    p2 = plist(idx2);
    [x2, y2, z2] = ind2sub(size(touchFace2), p2);
    coms(i, :) = mean([[x1, y1, z1]; [x2, y2, z2]], 1);
end

% to global coordinates
coms = bsxfun(@plus, coms, bboxSmall1(:,1)' - 1);
coms = round(coms);


end

