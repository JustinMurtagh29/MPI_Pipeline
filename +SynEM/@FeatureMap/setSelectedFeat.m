function setSelectedFeat(obj, selectedFeat)
%SETSELECTEDFEAT Determine the feature to be calculated if only a subset of
% features is required. This is only relevant for the 'direction' mode
% since in this case features always need to be calculated for both
% subsegments.
% INPUT selectedFeat: [Nx1] logical
%           Array of length obj.numFeatures selecting the features to
%           calculate.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isrow(selectedFeat)
    selectedFeat = selectedFeat';
end
if ~islogical(selectedFeat)
    error('SelectedFeat must be a logical array');
end
if length(selectedFeat) ~= obj.numFeatures
    error('selectedFeat must be specified for all features');
end

%reshape selected features
numSumStats = length(obj.quantiles) + sum(obj.moments);
numVols = 1 + 2.*obj.numSubvolumes;
nT = length(obj.featTexture);
featTextChannels = cellfun(@(x)x.numChannels, obj.featTexture);
featShapeNumFeat = cellfun(@(x)x.numFeatures, obj.featShape);

%transform selectedFeat to cell for each feature
selectedFeat = mat2cell(selectedFeat, ...
    [featTextChannels.*(numSumStats*numVols); featShapeNumFeat]);

%reshape texture features
selectedFeat(1:nT) = ...
    cellfun(@(x)reshape(x,numSumStats,numVols,[]), ...
    selectedFeat(1:nT),'UniformOutput',false);

%determine feature importance from input
switch obj.mode
    case 'single'
        %nothing to do
    case 'direction'
        %make sure that features are calculated for both both
        %subsegments
        subSegIntIdx = [1 reshape([3:2:2*(obj.numSubvolumes+1); ...
            2:2:2*(obj.numSubvolumes)],1,[])];
        selectedFeat(1:nT) = cellfun(@(x)x | x(:,subSegIntIdx,:), ...
            selectedFeat(1:nT), 'UniformOutput',false);
        if ~isempty(obj.featShapeDirDep)
            selectedFeat(nT+1:end) = cellfun(@(x, y)x | x(y), ...
                selectedFeat(nT+1:end), obj.featShapeDirDep, ...
                'UniformOutput', false);
        end
end

%save output in feature map
obj.selectedFeat = selectedFeat;
obj.numFeaturesSelected = sum(cellfun(@(x)sum(x(:)),selectedFeat));
end