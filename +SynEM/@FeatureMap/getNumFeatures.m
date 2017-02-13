function n = getNumFeatures( obj )
%GETNUMFEATURES Get the current number of selected features.
% OUTPUT n: int
%           Total number of features in the current feature map. This
%           corresponds to fm.numFeatures if no feature selection was done
%           or fm.numFeaturesSelected otherwise.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~isempty(obj.selectedFeat)
    n = obj.numFeaturesSelected;
else
    n = obj.numFeatures;
end


end

