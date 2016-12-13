function s = calculateSize(obj)
%CALCULATESIZE Caculate the numFeatures property, i.e. the total size of
% the feature map without feature selection.
% OUTPUT s: The size for the specified mode.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

numSumStats = length(obj.quantiles) + sum(obj.moments);
numVols = 1 + 2.*obj.numSubvolumes;
featTextChannels = sum(cellfun(@(x)x.numChannels, obj.featTexture));
featShape = sum(cellfun(@(x)x.numFeatures, obj.featShape));
obj.numFeatures = numSumStats*numVols*featTextChannels + featShape;
s = obj.numFeatures;

end