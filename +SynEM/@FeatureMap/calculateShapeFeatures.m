function X = calculateShapeFeatures(obj, i, coordVols)
%CALCULATESHAPEFEATURES Wrapper for shape feature calculation.
% INPUT i: Index of the shape feature.
%       coordVols: [Nx3] cell array containing the coordinates of the
%           interface surface and largest subvolumes as [Mx3] arrays.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

j = i + length(obj.featTexture);
selFeat = obj.selectedFeat{j}';
X = obj.featShape{i}.calculate( ...
    coordVols(:,obj.featShapeSubvols{i}),selFeat);

end
