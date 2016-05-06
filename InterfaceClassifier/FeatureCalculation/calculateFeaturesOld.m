function X = calculateFeaturesOld( raw, interfaces, featureMap )
%CALCULATEFEATURESOLD Calculate Features on output of interface calculation.
% INPUT raw: 3d raw data.
%       interfaces: see calculateInterfaces
%       featureMap: A previously created and maybe updated (feature
%                 selection) feature map.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%apply area threshold and load interfaces and raw data (compatibility)
indices = cellfun(@(x)length(x) > featureMap.areaThreshold,interfaces.surface);
interfaces.surface = interfaces.surface(indices);
interfaces.subseg = cellfun(@(x)x(indices,:),interfaces.subseg,'UniformOutput',false);
interfaceSurfaceList = interfaces.surface;
subsegmentsList = interfaces.subseg;
subsegmentsList = cat(2,subsegmentsList{1:length(subsegmentsList)});
raw = single(raw);

X = zeros(2*length(interfaceSurfaceList),sum(featureMap.numSelectedFeatures),'single');

fprintf('%s INFO - calculating features\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));
iter = 1;
for nFeature = 1:length(featureMap.type)
    if any(subsref(cell2mat(featureMap.featureIndices{nFeature}),struct('type','()','subs',{{':'}})))
        X(:,iter:iter + featureMap.numSelectedFeatures(nFeature) - 1) = calculateFeatureOld(featureMap.names{nFeature},featureMap.parameters{nFeature},raw,interfaceSurfaceList,subsegmentsList,featureMap.featureIndices{nFeature},featureMap.quantiles,featureMap.summaryStatistics,featureMap.border);
    end
    iter = iter + featureMap.numSelectedFeatures(nFeature);
end

end
