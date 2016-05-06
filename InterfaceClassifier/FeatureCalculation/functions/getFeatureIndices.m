function [featureMap, keptFeatureIndices] = getFeatureIndices( featureMap,varargin )
%GETFEATUREINDICES Determines the necessary features and updates the
%features map (necessary in this context means that if a feature is
%calculated for subsegments one, then it will also be calculated for
%subsegment two due to the need to interchange the two subsegments).
% varargin: Logical array which features to use (intended use: use imp>0
%           as featureIndices, where imp is the predictor importance
%           returned by a matlab ensemble)
% OUTPUT featureMap: The updated feature map.
%        keptFeatureIndices: Logical array for converting the original
%           feature matrix to the new feature set via
%           X(:,keptFeatureIndices)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isempty(varargin)
    featureIndices = true(sum(featureMap.numFeatures),1);
else
    featureIndices = varargin{1};
end

featureIndices = mat2cell(logical(featureIndices),featureMap.numFeatures);
featureSelectionIndices = cell(length(featureMap.type),1);
featureMap.numSelectedFeatures = zeros(length(featureMap.type),1);
numSubVolumes = length(featureMap.rinclude);
numPoolStat = length(featureMap.quantiles) + length(featureMap.summaryStatistics);
for k = 1:length(featureMap.type)
    switch featureMap.type{k}
        case 'texture'
            if featureMap.multiplier(k) > 1
                featureIndices{k} = reshape(featureIndices{k},numPoolStat,1+2*numSubVolumes,featureMap.multiplier(k));
                fIndices = mat2cell(featureIndices{k},size(featureIndices{k},1),size(featureIndices{k},2),ones(size(featureIndices{k},3),1));
                featureSelectionIndices{k} = fIndices(:);
            else
                featureSelectionIndices{k} = {reshape(featureIndices{k},numPoolStat,1+2*numSubVolumes)};
            end
            featureSelectionIndices{k} = cellfun(@(x) x | x(:,[1 reshape([3:2:2*(numSubVolumes+1);2:2:2*(numSubVolumes)],1,2*(numSubVolumes))]),featureSelectionIndices{k},'UniformOutput',false);
            
        case 'shape'
            if length(featureIndices{k}) == 3
                featureSelectionIndices{k} = {featureIndices{k} | featureIndices{k}([1 3 2])};
            else
                featureSelectionIndices{k} = featureIndices(k);
            end
    end
    featureMap.numSelectedFeatures(k) = sum(sum(sum(cell2mat(featureSelectionIndices{k}))));
end
featureMap.featureIndices = featureSelectionIndices;
keptFeatureIndices = cell2mat(cellfun(@(y)cell2mat(cellfun(@(x)reshape(x,numel(x),1),y,'UniformOutput',false)),featureMap.featureIndices,'UniformOutput',false));
featureMap.selectedFeatureInfo = featureMap.featureInfo(keptFeatureIndices);
featureMap.border = calculateFMBorder(featureMap.names,featureMap.parameters,featureSelectionIndices);

end