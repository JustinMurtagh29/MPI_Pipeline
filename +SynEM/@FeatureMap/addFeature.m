function addFeature(obj, feat, varargin)
%ADDFEATURE Adds a feature to the feature map.
%INPUT feat: An object from the Feature package.
%      varargin: Additional input arguments.
%           Shape Features:
%           varargin{1}: Integer vector containing the subvolumes
%               which are passed to the shape feature in the desired order.
%               1: interface surface
%               2,3: First and second subsegment of the largest subvolume.
%               e.g. [1 2 3] will pass the interface surface and the both
%                   subsegments of the largest subvolume as the vols
%                   argument to a shape feature.
%           varargin{2}: (Optional) Integer vector of 
%               length feat.numFeatures specifying how the features of the
%               feature vector need to be interchanged when the subvolumes
%               of an interface are interchanged, i.e. if X is the feature 
%               vector for one direction, then X(:,featShapeDirDep{i})
%               should be the  feature vector in the other direction.
%               This is only relevant for the 'direction' mode.
%               (Default: 1:feat.numFeatures)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

supClass = superclasses(feat);
switch supClass{1}
    case 'SynEM.Feature.TextureFeature'
        obj.featTexture{end+1,1} = feat;
        obj.border = max(obj.border, feat.border);
    case 'SynEM.Feature.ShapeFeature'
        obj.featShape{end+1,1} = feat;
        obj.featShapeSubvols{end+1,1} = varargin{1};
        if length(varargin) > 1
            obj.featShapeDirDep{end + 1,1} = varargin{2};
        else
            obj.featShapeDirDep{end + 1,1} = ...
                1:length(feat.numFeatures);
        end
    otherwise
        error(['Feature must be inherited from', ...
            'Feature.TextureFeature or Feature.ShapeFeature']);
end
calculateSize(obj);
end