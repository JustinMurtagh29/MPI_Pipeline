function arbitraryResliceForErrors(raw, seg, edge, colors, errorIdx, prefix)

% % Not that important, decide which segments to color in video & make sure
% % start and end segment always have the same color by always using setting
% % start segment as object 1 and end object as object 2
segNew = zeros(size(seg));
segNew(seg == edge(1)) = 1;
segNew(seg == edge(2)) = 2;
props = regionprops(imdilate(segNew == 1, ones(3,3,3)) & imdilate(segNew == 2, ones(3,3,3)));
[~,maxIdx] = max([props(:).Area]);
props = props(maxIdx);
errorCenter = round(props.Centroid([2 1 3]));
props1 = regionprops(segNew == 1);
props2 = regionprops(segNew == 2);
direction = round(props2.Centroid([2 1 3]) - props1.Centroid([2 1 3]));
direction = direction ./ norm(direction, 2);
[~, orth1, orth2] = findRandomOrthogonal(direction);
colors = colors(1:max(segNew(:)),:);

% set/extract some needed parameters
frameSize = [200 200]; % same as in mission.json
renderedFrames = [40 40]; % setting from levelcreator (how many frames before and after problem location)
voxelSize = [11.24 11.24 28]; % from settings.json (anisotropy)
voxelSizeAfter = [11.24 11.24 11.24]; % size of voxel after interpolation (in video)

% Find out how much data we need around error center (calculations done
% on small cube for speed purposes)
border = sqrt(2) .* max(voxelSizeAfter .* [frameSize/2 max(renderedFrames)]); % border needed in nm
% sqrt(2) ue to maximal effect if reslice is 45 degress wrt original
% orientation in plane
border = ceil(repmat(border,1,3) ./ voxelSize); % divide by voxel size of original data and round up

% Calculate bounding box from border
bbox = [errorCenter - border; errorCenter + border]; %bounding box in dataset coordinates

% cut local stack for interpolation & rendering according to bbox
rawLocal = raw(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));
segLocal = segNew(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));

% Define grid for original values (rawLocal & segLocal & mitoLocal),
% which has size 2* border+1 and is centered on the error center
[X,Y,Z] = meshgrid(-border(1):border(1),-border(2):border(2),-border(3):border(3)); % grid in dataset coordinates [voxel]
X = X .* voxelSize(1); % grid in dataset coordinates [nm]
Y = Y .* voxelSize(2);
Z = Z .* voxelSize(3);
% Construct rotation matrix from orthogonal basis (one of which is direction)
A = [orth1; orth2; direction];
% Define grid for interpolating values (size according to settings from
% viewport & levelcreator
[XI, YI, ZI] = meshgrid(-frameSize(1)/2:frameSize(1)/2, -frameSize(2)/2:frameSize(2)/2, -renderedFrames(1):renderedFrames(2)); %new grid [voxel]
% Linearize representation of points and scale according to desired
% voxel size
coords = bsxfun(@times, [XI(:) YI(:) ZI(:)], voxelSizeAfter); % new grid transform to [nm] including making 3D Matrices X,Y,Z to 1D for rotation
coords = coords*A; % rotate new grid according to ONB
XI = reshape(coords(:,1),size(XI)); % Undo linearization
YI = reshape(coords(:,2),size(YI));
ZI = reshape(coords(:,3),size(ZI));
% Do the interpolation from original grid (X,Y,Z) indicating location
% of (raw/seg/mito)Local values onto rotated grid with different
% anisotropy (all grids are in nm)
rawLocalNewDirection = interp3(X, Y, Z, rawLocal, XI, YI, ZI, 'cubic');
segLocalNewDirection = interp3(X, Y, Z, single(segLocal), XI, YI, ZI, 'nearest');
% Make movies of stacks (no more processing done in there, just
% displaying frame by frame and capturing)
makeTaskMovieSetColors(uint16(segLocalNewDirection), rawLocalNewDirection, ['C:\Users\mberning\Desktop\classifier\errors\' prefix num2str(errorIdx, '%.5i') '.avi'], colors);

function [randomVec, orth1, orth2] = findRandomOrthogonal(v)
% returns random orthogonal vector and and the 2nd and 3rd vector for ONB
v = v ./ norm(v);
if all(v == [0 0 1])
    orth1 = [1 0 -1].*v([3 2 1]);
else
    orth1 = [1 -1 0].*v([2 1 3]);
end
orth1 = orth1 ./ norm(orth1);
orth2 = cross(v, orth1);
mix = 2 * (rand(2,1) -0.5);
randomVec = mix(1)*orth1+mix(2)*orth2;
randomVec = randomVec ./ norm(randomVec);
end

end
