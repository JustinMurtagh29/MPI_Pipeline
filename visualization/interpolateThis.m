function [rawForLevelcreator, segForLevelcreator] = interpolateThis(raw, seg, errorCenter, direction, frameSize, renderedFrames, voxelSizeBefore, voxelSizeAfter)
% Create new stack which is an arbitrary reslice of the old stack
% Inputs:
% 1) 3D matrix with raw values (cut out of raw data)
% 2) 3D matrix with correspoding segmentation (layer/segmentation)
% 3) 3D vector: error center (missions.json)
% 4) 3D vector: direction of flight (missions.json)
% 4) 2D vector: requested size of output frame in voxel (set in levelcreator)
% 5) 2D vector: frames before/after error center in voxel (set in levelcreator)
% 7) Voxel size raw data (see settings.json -> voxelSize)
% 8) Voxel size in output stack (that is the number I was interested in last week)

%% 1st: Cut out a substack
% Calculate how much context around error center we need in nm
% sqrt(2) comes from 45 degree rotation (pythagoras)
% /2 comes from the fact that we are calculating border on each size of error center
% inner max is needed if not centered on error center in time direction of the stack
% outer max to take into account all possible rotations
border = sqrt(2) .* max(voxelSizeAfter .* [frameSize/2 max(renderedFrames)]);
% now going back to voxel space for cutting (we obviously need to cut less voxel in z-diretcion
border = ceil(repmat(border,1,3) ./ voxelSizeBefore);
% Bounding Box to cut out of global dataset (still everything in voxel)
bbox = [errorCenter - border; errorCenter + border];
% cut local stacks for interpolation & rendering 
rawLocal = raw(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));
segLocal = seg(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));

%% 2nd: Transfer voxel measures in JSON to nm space
% calculate direction in nm
direction = direction .* voxelSizeBefore; % THIS IS ONLY NECESSARY BECAUSE JSON direction is still in voxel
% Normalize direction again
direction = direction ./ norm(direction);
% calculate position of error center in nm
errorCenter = errorCenter .* voxelSizeBefore;

%% 3rd: define nm-grid of the stack rawLocal & segLocal
% Define grid for raw values (where are the raw grayvalues/segmentation in nm space centered on error center)
[X,Y,Z] = meshgrid(-border(1):border(1),-border(2):border(2),-border(3):border(3));
% Scale voxel grid to nm & translate so that it fits on error center
X = (X .* voxelSizeBefore(1) + errorCenter(1));
Y = (Y .* voxelSizeBefore(2) + errorCenter(2));
Z = (Z .* voxelSizeBefore(3) + errorCenter(3));

%% 4th: Construct rotation matrix
% Construct orthogonal basis (last of which is direction, so that direction will be z-direction in new stack)
orth1 = [-1 1 0].*direction([2 1 3]);
orth1 = orth1 ./ norm(orth1);
orth2 = -cross(direction , orth1);
% Construct rotation matrix from orthogonal basis
A = [orth1; orth2; direction];

%% 5th: What points in nm-space do we want to know the interpolated grey values of
% Define grid for interpolating values
[XI, YI, ZI] = meshgrid(-frameSize(1)/2:frameSize(1)/2, -frameSize(2)/2:frameSize(2)/2, -renderedFrames(1):renderedFrames(2));
coords = [XI(:) YI(:) ZI(:)];
coords = coords .* repmat(voxelSizeAfter, size(coords, 1),1);
% Rotation of new grid
coords = coords/A;
% Translation to error center
coords = coords + repmat(errorCenter, size(coords,1),1);
% Split and reshape arrays to fit interp3
XI = reshape(coords(:,1),size(XI));
YI = reshape(coords(:,2),size(YI));
ZI = reshape(coords(:,3),size(ZI));

%% 6th: interpolate
% Do the interpolation from original onto new (rotated & translated) grid
rawForLevelcreator = interp3(X, Y, Z, rawLocal, XI, YI, ZI, 'cubic');
segForLevelcreator = interp3(X, Y, Z, single(segLocal), XI, YI, ZI, 'nearest');

end
