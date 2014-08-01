function result = calculateLocalCorrespondences(cubeCoords1, cubeCoords2, segFile1, segFile2, localBorder, nrPlanesToCompare, saveFolder )
% calculateLocalCorrespondences -> Find correspondence between local segmentations by comparing overlap around face between cubes

% Load both segmentations
cube1 = load(segFile1);
cube2 = load(segFile2);

% How the cubes are arranged wrt to each other
coordShift = cubeCoords2' - cubeCoords1';

% some needed sizes?
outerCubeSize = size(cube1.seg)'; % this function needs constant sized bounding boxes to work properly
innerCubeSize = outerCubeSize + localBorder(:,1) - localBorder(:,2);
bbox = cat(2,[1;1;1], outerCubeSize) - localBorder;

% Calculating bounding boxes for segmentations to compare against (@Alex: i think the 'any' in your code was an error?)
% Should be a 'cleaner' / more general approach to do this?
bbox1 = bbox; % inner bounding box for cube1
bbox1(:,2) = bbox(:,2) - coordShift .* innerCubeSize; % align with one face
bbox1(:,1) = bbox1(:,1) - nrPlanesToCompare .* coordShift .* [1;1;1]; % expand along dimension of inter-cube displacement to create overlap
bbox1(:,2) = bbox1(:,2) + nrPlanesToCompare .* coordShift .* [1;1;1]; % same for 2nd dimension
bbox1 = bbox1 + repmat((coordShift > 0) .* innerCubeSize, [1 2]); % shift to 'upward' face in case coordShift is positive
% same for 2nd cube (slightly altered last line due to asymetric coordShift wrt cubeCoords)
bbox2 = bbox; % inner bounding box for cube2
bbox2(:,2) = bbox(:,2) - coordShift .* innerCubeSize; % align with one face
bbox2(:,1) = bbox2(:,1) - nrPlanesToCompare .* coordShift .* [1;1;1]; % expand along dimension of inter-cube displacement to create overlap
bbox2(:,2) = bbox2(:,2) + nrPlanesToCompare .* coordShift .* [1;1;1]; % same for 2nd dimension
bbox2 = bbox2 + repmat((coordShift < 0) .* innerCubeSize, [1 2]); % shift to 'upward' face in case coordShift is negative

% Get segmentations to compare for correspondence finding
overlap1 = cube1.seg(bbox1(1,1):bbox1(1,2), bbox1(2,1):bbox1(2,2), bbox1(3,1):bbox1(3,2));
overlap2 = cube2.seg(bbox2(1,1):bbox2(1,2), bbox2(2,1):bbox2(2,2), bbox2(3,1):bbox2(3,2));
clear cube1 cube2;

% create one dimensional array to compare both overlaps
correspondences = [reshape(overlap1, [numel(overlap1) 1 1]) reshape(overlap2, [numel(overlap2) 1 1])];
% remove correspondences to 0 (background)
correspondences(any(correspondences == 0,2),:) = [];

% Make dimension orthogonal to plane last dimension in array
permuteVector = [1 2 3];
permuteVector(permuteVector == find(coordShift)) = 3;
permuteVector(3) = find(coordShift);
overlap1 = permute(overlap1, permuteVector);
overlap2 = permute(overlap2, permuteVector);

% Find all segments that appear in z=2 of overlap1 but not z=1 and all segments that appear z=1 in overlap2 but not z=2
o1z1 = unique(reshape(overlap1(:,:,1), [numel(overlap1)/2 1 1]));
o1z2 = unique(reshape(overlap1(:,:,2), [numel(overlap1)/2 1 1]));
o1ToRemove = setdiff(o1z2,o1z1);
o2z1 = unique(reshape(overlap2(:,:,1), [numel(overlap2)/2 1 1]));
o2z2 = unique(reshape(overlap2(:,:,2), [numel(overlap2)/2 1 1]));
o2ToRemove = setdiff(o2z1,o2z2);
% Remove those
for i=1:length(o1ToRemove)
    correspondences(correspondences(:,1) == o1ToRemove(i),:) = [];
end
for i=1:length(o2ToRemove)
    correspondences(correspondences(:,2) == o2ToRemove(i),:) = [];
end

% find unique correspondences
uniqueCorrespondences = unique(correspondences, 'rows');

% create output struct
result.cubeCoords1 = cubeCoords1;
result.cubeCoords2 = cubeCoords2;
result.overlap1 = overlap1;
result.overlap2 = overlap2;
result.correspondences = uniqueCorrespondences;

% Create main correspondence folder if necessary
if ~exist(saveFolder, 'dir')
	mkdir(saveFolder);
end

save([saveFolder num2str(cubeCoords1, '%.2i') num2str(cubeCoords2, '%.2i') '.mat'], 'result');

end

