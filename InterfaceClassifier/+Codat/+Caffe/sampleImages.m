function [ X, y, coords ] = sampleImages( raw, target, imSize, numClassInstances, augmentData )
%SAMPLEIMAGES Sample images from raw.
% INPUT raw: 3D raw data stack.
%       target: 3D binary target stack which centered in raw. The size
%               difference of raw and target must be enough such that one
%               can fit the imSize around every target voxel in raw.
%       imSize: Size of the resulting images for each target voxel (size
%               must be odd in each dimension.
%       numClassInstances: How many images to sample per class.
%       augmentData: Flag specifying whether to use data augmentation on
%               the images (currently 4x larger dataset).
% OUTPUT X: Image array of shape length x width x instances.
%        y: Labels for corresponding instances in X. True labels correspond
%           to target == 1 labels.
%        coords: Center of corresponding instance images in X as linear
%           indices in raw.

border = (size(raw) - size(target))/2;
if any((imSize - 1)/2 > [border(1), border(2)])
    error('Image size too larger for border');
end
imHSize = (imSize - 1)/2;
posLabels = find(target == 1);
negLabels = find(target == 0);
    
if numClassInstances > length(posLabels)
    warning('Not enough labels for requested number of class instances. Reducing to maximal possible number');
    numClassInstances = length(posLabels);
end

X = zeros([imSize,2*numClassInstances],'like',raw);
y = false(1,2*numClassInstances);
coords = zeros(1,2*numClassInstances);
posIdx = randperm(length(posLabels));
negIdx = randperm(length(negLabels));
for i = 1:numClassInstances
    %pos instances
    [u,v,w] = ind2sub(size(target),posLabels(posIdx(i)));
    u = u + border(1); v = v + border(2); w = w + border(3);
    X(:,:,2*i -1) = raw(u - imHSize(1):u + imHSize(1),v - imHSize(2):v + imHSize(2),w);
    y(2*i - 1) = true;
    coords(2*i - 1) = posLabels(posIdx(i));
    %neg instances
    [u,v,w] = ind2sub(size(target),negLabels(negIdx(i)));
    u = u + border(1); v = v + border(2); w = w + border(3);
    X(:,:,2*i) = raw(u - imHSize(1):u + imHSize(1),v - imHSize(2):v + imHSize(2),w);
    y(2*i) = false;
    coords(2*i) = negLabels(negIdx(i));
end

if augmentData
    original = X;
    for k = 1:3
        X = cat(3,X,rot90(original,k));
    end
    y = repmat(y,1,4);
    coords = repmat(coords,1,4);
end


end

