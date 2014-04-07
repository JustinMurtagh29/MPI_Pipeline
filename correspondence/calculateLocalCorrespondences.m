function result = calculateLocalCorrespondences(cubeCoords1, cubeCoords2, segFile1, segFile2, localBorder, nrPlanesToCompareOnEachSide )
% calculateLocalCorrespondences -> Find correspondence between local segmentations by comparing overlap around face between cubes

% How the cubes are arranged wrt to each other
coordShift = cubeCoords2 - cubeCoords1;

% some needed sizes?
outerCubeSize = size(cube1.seg); % this function needs constant sized bounding boxes to work properly
innerCubeSize = outerCubeSize + localBorder(:,1) - localBorder(:,2);

% Calculating bounding boxes for segmentations to compare against (@Alex: i think the 'any' in your code was an error?)
shiftWithinInnerCube =  (innerCubeSize - nrPlanesToCompareOnEachSize).*coordShift;
bbox1(:,1) = abs(localBorder(:,1)) + [1 1; 1 1; 1 1] + ;
bbox1(:,1) = outerCubeSize - localBorder(:,2) + (innerCubeSize - nrPlanesToCompareOnEachSide).*coordShift
bbox1(i,1) = -localBorder(:,1) + (innerCubeSize.*coordS- nrPlanesToCompareOnEachSide).*coordShift; 
bbox1(i,2) = outerCubeSizelocalBorder(:,2) + innerCubeSize.*coordShift - 
cubeSize(i)*3/4 + any(coordShift(i));
bbox2(i,1) = cubeSize(i)/4 - any(coordShift(i));
bbox2(i,2) = cubeSize(i)*3/4 - any(coordShift(i))*0.5*cubeSize(i);

% Load both segmentations
cube1 = load(segFile1);
cube2 = load(segFile2);

% Get segmentations to compare for correspondence finding
overlap1 = cube1(bbox1(1,1):bbox1(1,2), bbox1(2,1):bbox1(2,2), bbox1(3,1):bbox1(3,2));
overlap2 = cube2(bbox2(1,1):bbox2(1,2), bbox2(2,1):bbox2(2,2), bbox2(3,1):bbox2(3,2));

% create one dimensional array to compare both overlaps
overlap1 = reshape(overlap1, [numel(overlap1) 1 1]);
overlap2 = reshape(overlap2, [numel(overlap2) 1 1]);

%create output struct
ind = [ind1 ind2];
[ind(:,1), idx] = sort(ind(:,1));
ind(:,2) = ind(idx,2);
ind(any(ind == 0,2),:) = [];
result.ind = uint32(unique(ind, 'rows'));
result.ind(:,1) = result.ind(:,1) + 10000000 * cubeCoords1(1) +  1000000 * cubeCoords1(2) +  100000 * cubeCoords1(3); 
result.ind(:,2) = result.ind(:,2) + 10000000 * cubeCoords2(1) +  1000000 * cubeCoords2(2) +  100000 * cubeCoords2(3);

result.cubeCoords1 = cubeCoords1;
result.cubeCoords2 = cubeCoords2;

result.overlap1 = overlap1;
result.overlap2 = overlap2;

end
