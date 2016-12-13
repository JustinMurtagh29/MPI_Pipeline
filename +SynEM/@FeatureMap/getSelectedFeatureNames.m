function names = getSelectedFeatureNames(obj)
%GETSELECTEDFEATURENAMES Get the names for the currently selected features.
% OUTPUT names: [Nx1] cell array of length obj.numFeaturesSelected
%           containing the names of the currently selected features.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isempty(obj.names)
    obj.createNameStrings();
end
if isempty(obj.selectedFeat)
    names = obj.names;
else
    names = obj.names(cell2mat(cellfun( ...
        @(x)x(:),obj.selectedFeat, 'UniformOutput', false)));
end
end